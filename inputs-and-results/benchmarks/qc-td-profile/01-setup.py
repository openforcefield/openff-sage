import copy
import functools
import json
from collections import defaultdict
from multiprocessing import Pool
from typing import Dict, Tuple

import click
from openff.qcsubmit.results import TorsionDriveResultCollection
from openff.qcsubmit.results.filters import RecordStatusFilter
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from openmmforcefields.generators import GAFFTemplateGenerator
from qcportal import FractalClient
from qcportal.models import TorsionDriveRecord
from qcportal.models.records import RecordStatusEnum
from simtk import unit
from simtk.openmm import openmm
from simtk.openmm.app import ForceField as OMMForceField
from tqdm import tqdm


@functools.lru_cache()
def _parameterize_molecule(
    mapped_smiles: str,
    force_field_path: str,
    force_field_type: str,
) -> openmm.System:
    """Parameterize a particular molecules with a specified force field."""

    off_molecule = Molecule.from_mapped_smiles(mapped_smiles)

    if force_field_type.lower() == "smirnoff":
        force_field = ForceField(force_field_path, allow_cosmetic_attributes=True)

        if "Constraints" in force_field.registered_parameter_handlers:
            force_field.deregister_parameter_handler("Constraints")

        return force_field.create_openmm_system(off_molecule.to_topology())

    elif force_field_type.lower() == "gaff":

        force_field = OMMForceField()

        generator = GAFFTemplateGenerator(
            molecules=off_molecule, forcefield=force_field_path
        )

        force_field.registerTemplateGenerator(generator.generator)

        return force_field.createSystem(
            off_molecule.to_topology().to_openmm(),
            nonbondedCutoff=0.9 * unit.nanometer,
            constraints=None,
        )

    raise NotImplementedError(
        "Only SMIRNOFF and GAFF force fields are currently supported."
    )


def _evaluate_energy(openmm_system: openmm.System, coordinates: unit.Quantity) -> float:
    """Returns the evaluated energy of the conformer in kcal / mol."""

    integrator = openmm.VerletIntegrator(0.001 * unit.femtoseconds)

    platform = openmm.Platform.getPlatformByName("Reference")
    openmm_context = openmm.Context(openmm_system, integrator, platform)

    openmm_context.setPositions(coordinates.value_in_unit(unit.nanometers))
    state = openmm_context.getState(getEnergy=True)

    potential_energy = state.getPotentialEnergy()

    return potential_energy.value_in_unit(unit.kilocalories_per_mole)


def _minimise_structure(
    openmm_system: openmm.System,
    coordinates: unit.Quantity,
    fixed_indices: Tuple[int, ...],
) -> unit.Quantity:
    """Minimizes a set of coordinates using a specified potential energy function."""

    openmm_system = copy.deepcopy(openmm_system)

    # Constrain the fixed atoms.
    for index in fixed_indices:
        openmm_system.setParticleMass(index, 0.0)

    integrator = openmm.VerletIntegrator(0.001 * unit.femtoseconds)

    platform = openmm.Platform.getPlatformByName("Reference")

    openmm_context = openmm.Context(openmm_system, integrator, platform)
    openmm_context.setPositions(coordinates.value_in_unit(unit.nanometers))

    openmm.LocalEnergyMinimizer.minimize(openmm_context)

    state: openmm.State = openmm_context.getState(getPositions=True)
    minimised_coordinates = state.getPositions()

    return minimised_coordinates


def _compute_profile_and_residual(
    force_field_path: str,
    force_field_type: str,
    qc_record: TorsionDriveRecord,
    molecule: Molecule,
) -> Dict[Tuple[int, ...], Tuple[float, float, float]]:
    """Convert a QC record and its associated molecule into an OE molecule which
    has been tagged with the associated SMILES, final energy and record id."""

    grid_energies = qc_record.get_final_energies()

    grid_conformers = {
        tuple(json.loads(grid_id)): conformer
        for grid_id, conformer in zip(
            molecule.properties["grid_ids"], molecule.conformers
        )
    }

    grid_ids = sorted(grid_conformers, key=lambda x: x[0])
    dihedral_indices = qc_record.keywords.dihedrals[0]

    openmm_system = _parameterize_molecule(
        molecule.to_smiles(isomeric=True, mapped=True),
        force_field_path,
        force_field_type,
    )
    residual_system = copy.deepcopy(openmm_system)

    # Zero out the contribution of the driven torsion.
    torsion_force = [
        force
        for force in residual_system.getForces()
        if isinstance(force, openmm.PeriodicTorsionForce)
    ][0]

    for torsion_index in range(torsion_force.getNumTorsions()):

        i, j, k, l, periodicity, phase, _ = torsion_force.getTorsionParameters(
            torsion_index
        )

        if sorted([i, j]) != sorted(dihedral_indices[1:3]):
            continue

        torsion_force.setTorsionParameters(
            torsion_index, i, j, k, l, periodicity, phase, 0.0
        )

    # Evaluate the energy for each grid id.
    energies: Dict[Tuple[int, ...], Tuple[float, float, float]] = dict()

    lowest_qm_energy = None
    lowest_qm_energy_grid_id = None

    lowest_mm_residual = None

    for grid_id in grid_ids:

        coordinates = _minimise_structure(
            openmm_system, grid_conformers[grid_id], dihedral_indices
        )

        qm_energy = (
            grid_energies[grid_id] * unit.hartree * unit.AVOGADRO_CONSTANT_NA
        ).value_in_unit(unit.kilocalories_per_mole)

        mm_energy = _evaluate_energy(openmm_system, coordinates)
        mm_residual = qm_energy - _evaluate_energy(residual_system, coordinates)

        if lowest_qm_energy is None or qm_energy < lowest_qm_energy:

            lowest_qm_energy = qm_energy
            lowest_qm_energy_grid_id = grid_id

        if lowest_mm_residual is None or mm_residual < lowest_mm_residual:
            lowest_mm_residual = mm_residual

        energies[grid_id] = (qm_energy, mm_energy, mm_residual)

    energies = {
        grid_id: (
            qm_energy - energies[lowest_qm_energy_grid_id][0],
            mm_energy - energies[lowest_qm_energy_grid_id][1],
            mm_residual - lowest_mm_residual,
        )
        for grid_id, (qm_energy, mm_energy, mm_residual) in energies.items()
    }

    return energies


@click.command()
@click.option(
    "-d",
    "--data-set",
    "data_set_name",
    default="OpenFF-benchmark-ligand-fragments-v1.0",
    show_default=True,
    type=click.STRING,
    help="The name of the QCArchive torsion drive collection to computed the "
    "QM, MM energies and residuals for.",
)
@click.option(
    "-ff",
    "--force-fields",
    "force_field_input_path",
    default="01-force-fields.json",
    show_default=True,
    type=click.Path(exists=True, dir_okay=False),
    help="The path to the JSON file which specifies which force fields to use to "
    "compute the MM torsion profile and residual.",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    default="01-residuals.json",
    show_default=True,
    type=click.Path(exists=False, dir_okay=False),
    help="The path to save the QM profile, MM profile and MM residuals to.",
)
def main(data_set_name, force_field_input_path, output_path):

    with open(force_field_input_path) as file:
        force_fields = json.load(file)

    print("1a) Parsing collection")

    client = FractalClient()

    result_collection = TorsionDriveResultCollection.from_server(
        client=client,
        datasets=data_set_name,
        spec_name="default",
    )

    print("1b) Filtering unfinished results")

    result_collection = result_collection.filter(
        RecordStatusFilter(status=RecordStatusEnum.complete)
    )

    print("1c) Retrieving QC records")

    records_and_molecules = result_collection.to_records()

    print("1d) Computing the residuals")

    profile_and_residuals = defaultdict(lambda: defaultdict(dict))
    record_id_to_smiles = {}

    with Pool(processes=len(force_fields)) as pool:

        for qc_record, off_molecule in tqdm(records_and_molecules):

            force_field_labels = [*force_fields]
            force_field_path_and_types = [
                (force_fields[label]["path"], force_fields[label]["type"])
                for label in force_field_labels
            ]

            force_field_energies = pool.starmap(
                functools.partial(
                    _compute_profile_and_residual,
                    qc_record=qc_record,
                    molecule=off_molecule,
                ),
                force_field_path_and_types,
            )

            for label, energies in zip(force_field_labels, force_field_energies):

                for grid_id, energy in energies.items():
                    profile_and_residuals[qc_record.id][grid_id[0]][label] = energy

            # Track which SMILES (with the driven dihedral tagged) corresponds to which
            # record id.
            dihedral_indices = qc_record.keywords.dihedrals[0]

            off_molecule_copy = copy.deepcopy(off_molecule)
            off_molecule_copy.properties["atom_map"] = {
                j: i + 1 for i, j in enumerate(dihedral_indices)
            }
            record_id_to_smiles[qc_record.id] = off_molecule_copy.to_smiles(mapped=True)

    print(
        f"1e) Writing energies and residuals to JSON (N={len(records_and_molecules)})"
    )

    with open(output_path, "w") as file:

        json.dump(
            {"record_smiles": record_id_to_smiles, "energies": profile_and_residuals},
            file,
        )


if __name__ == "__main__":
    main()
