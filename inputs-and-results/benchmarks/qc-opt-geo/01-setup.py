import logging
import warnings
from collections import defaultdict
from multiprocessing import get_context
from typing import Tuple

import click
import qcelemental
from openeye import oechem
from openff.qcsubmit.results import OptimizationResultCollection
from openff.qcsubmit.results.filters import (
    CMILESResultFilter,
    ConnectivityFilter,
    RecordStatusFilter,
    SMILESFilter,
)
from openff.toolkit.topology import Molecule
from openff.toolkit.topology.molecule import SmilesParsingError
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils import (
    GLOBAL_TOOLKIT_REGISTRY,
    OpenEyeToolkitWrapper,
    UndefinedStereochemistryError,
)
from qcportal import FractalClient
from qcportal.models.records import RecordStatusEnum
from tqdm import tqdm

N_PROCESSES = 16


class InvalidCMILESFilter(CMILESResultFilter):
    def _filter_function(self, result) -> bool:

        try:
            Molecule.from_mapped_smiles(result.cmiles, allow_undefined_stereo=True)
        except (ValueError, SmilesParsingError):
            return False

        return True


def _can_parameterize(smiles: str) -> Tuple[str, bool]:

    try:

        for toolkit in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:

            if isinstance(toolkit, OpenEyeToolkitWrapper):
                continue

            GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(toolkit)

        molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
        force_field = ForceField("openff-1.3.0.offxml")

        force_field.create_openmm_system(molecule.to_topology())

    except:
        return smiles, False

    return smiles, True


def _process_molecule(record_and_molecule) -> oechem.OEMol:
    """Convert a QC record and its associated molecule into an OE molecule which
    has been tagged with the associated SMILES, final energy and record id."""

    record, molecule = record_and_molecule

    oe_molecule = molecule.to_openeye()
    oechem.OE3DToInternalStereo(oe_molecule)

    final_energy = record.get_final_energy() * qcelemental.constants.hartree2kcalmol

    # add name and energy tag to the mol
    oechem.OESetSDData(oe_molecule, "SMILES QCArchive", molecule.to_smiles())
    oechem.OESetSDData(oe_molecule, "Energy QCArchive", str(final_energy))
    oechem.OESetSDData(oe_molecule, "Record QCArchive", str(record.id))

    return oe_molecule


@click.command()
@click.argument(
    "data_set",
    nargs=1,
    type=click.STRING,
)
def main(data_set):
    """Process and store optimized QC geometries from a QCArchive dataset.

    The DATA_SET should likely be one of:

    \b
    * "OpenFF Full Optimization Benchmark 1"
    * "OpenFF Industry Benchmark Season 1 v1.0"

    The processed molecules tagged with information from the QC record, including the
    CMILES and QC energy, will be stored in a new `01-processed-qm.sdf` file and
    additionally information about the included records will be stored in
    `01-processed-qm.json`.
    """

    warnings.filterwarnings("ignore")
    logging.getLogger("openff.toolkit").setLevel(logging.ERROR)

    # Make sure we consistently only use OE in this script
    for toolkit in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits:
        if isinstance(toolkit, OpenEyeToolkitWrapper):
            continue
        GLOBAL_TOOLKIT_REGISTRY.deregister_toolkit(toolkit)

    print("1a) Parsing collection")

    client = FractalClient()

    result_collection = OptimizationResultCollection.from_server(
        client=client,
        datasets=data_set,
        spec_name="default",
    )

    print("1b) Filtering unfinished / undesirable results")

    result_collection = result_collection.filter(
        InvalidCMILESFilter(),
        RecordStatusFilter(status=RecordStatusEnum.complete),
        ConnectivityFilter(),
    )

    _, molecules = zip(*result_collection.to_records())

    unique_molecules = set()

    for molecule in molecules:

        # Re-perceive the stereochemistry from the final conformer.
        oe_molecule = molecule.to_openeye()
        oechem.OE3DToInternalStereo(oe_molecule)

        try:
            molecule = Molecule.from_openeye(oe_molecule)
        except UndefinedStereochemistryError:
            print(f"skipping {molecule.to_smiles()} - un-perceivable stereo")
            continue

        unique_molecules.add(molecule.to_smiles(isomeric=True, mapped=False))

    with get_context("spawn").Pool(processes=N_PROCESSES) as pool:

        filtered_smiles = {
            smiles
            for smiles, should_retain in tqdm(
                pool.imap_unordered(_can_parameterize, unique_molecules),
                total=len(unique_molecules),
            )
            if should_retain
        }

    result_collection = result_collection.filter(
        SMILESFilter(smiles_to_include=[*filtered_smiles]),
    )

    print("1e) Grouping molecules")

    records_and_molecules = result_collection.to_records()

    grouped_molecules = defaultdict(list)

    for record, molecule in records_and_molecules:

        molecule = molecule.canonical_order_atoms()

        smiles = molecule.to_smiles(isomeric=False, explicit_hydrogens=True)
        grouped_molecules[smiles].append((record, molecule))

    print("1f) Processing molecules")

    processed_oe_molecules = [
        _process_molecule(record_and_molecule)
        for record_and_molecule in records_and_molecules
    ]

    print("1g) Writing processed molecules to sdf")

    output_steam = oechem.oemolostream("01-processed-qm.sdf")

    final_record_ids = set()

    for i, oe_molecule in enumerate(processed_oe_molecules):

        final_record_ids.add(oechem.OEGetSDData(oe_molecule, "Record QCArchive"))

        oe_molecule.SetTitle(f"full_{i + 1}")
        oechem.OEWriteConstMolecule(output_steam, oe_molecule)

    output_steam.close()

    result_collection.entries[client.address] = [
        entry
        for entry in result_collection.entries[client.address]
        if entry.record_id in final_record_ids
    ]

    with open("01-processed-qm.json", "w") as file:
        file.write(result_collection.json())


if __name__ == "__main__":
    main()
