import json
from typing import Any, Dict

import click
import numpy
import numpy as np
import openeye.oechem as oechem
import pandas
from geometric.internal import (
    Angle,
    Dihedral,
    Distance,
    OutOfPlane,
    PrimitiveInternalCoordinates,
)
from geometric.molecule import Molecule as GeometricMolecule
from openff.toolkit.topology import Molecule
from rdkit.Chem import TorsionFingerprints
from simtk import unit

# Define the RMSD calculation parameters
RMSD_AUTOMORPH = True  # take into acct symmetry related transformations
RMSD_HEAVY_ONLY = False  # do consider hydrogen atoms for automorphisms
RMSD_OVERLAY = True  # find the lowest possible RMSD

def compute_rmsd(ref, tar, v_periodic=None):
    """
    Compute the RMSD between two arrays, supporting periodic difference
    Code adapted from FORCEBALANCE
    Source: https://github.com/leeping/forcebalance/blob/master/src/opt_geo_target.py
    """
    assert len(ref) == len(tar), 'array length must match'
    n = len(ref)
    if n == 0: return 0.0
    if v_periodic is not None:
        diff = periodic_diff(ref, tar, v_periodic)
    else:
        diff = ref - tar
    rmsd = np.sqrt(np.sum(diff**2) / n)
    return rmsd

def periodic_diff(a, b, v_periodic):
    ''' convenient function for computing the minimum difference in periodic coordinates
    Code adapted from FORCEBALANCE
    Source: https://github.com/leeping/forcebalance/blob/master/src/opt_geo_target.py

    Parameters
    ----------
    a: np.ndarray or float
        Reference values in a numpy array
    b: np.ndarray or float
        Target values in a numpy arrary
    v_periodic: float > 0
        Value of the periodic boundary
    Returns
    -------
    diff: np.ndarray
        The array of same shape containing the difference between a and b
        All return values are in range [-v_periodic/2, v_periodic/2),
        "( )" means exclusive, "[ ]" means inclusive
    Examples
    -------
    periodic_diff(0.0, 2.1, 2.0) => -0.1
    periodic_diff(0.0, 1.9, 2.0) => 0.1
    periodic_diff(0.0, 1.0, 2.0) => -1.0
    periodic_diff(1.0, 0.0, 2.0) => -1.0
    periodic_diff(1.0, 0.1, 2.0) => -0.9
    periodic_diff(1.0, 10.1, 2.0) => 0.9
    periodic_diff(1.0, 9.9, 2.0) => -0.9
    '''
    assert v_periodic > 0
    h = 0.5 * v_periodic
    return (a - b + h) % v_periodic - h


def _compute_internal_coordinate_rmsd(
    molecule: Molecule,
    qm_conformer: unit.Quantity,
    mm_conformer: unit.Quantity,
) -> Dict[str, float]:

    qm_conformer = qm_conformer.value_in_unit(unit.angstrom)
    mm_conformer = mm_conformer.value_in_unit(unit.angstrom)

    geo_molecule = GeometricMolecule()
    geo_molecule.Data = {
        "resname": ["UNK"] * molecule.n_atoms,
        "resid": [0] * molecule.n_atoms,
        "elem": [atom.element.symbol for atom in molecule.atoms],
        "bonds": [(bond.atom1_index, bond.atom2_index) for bond in molecule.bonds],
        "name": molecule.name,
        "xyzs": [qm_conformer],
    }

    internal_coordinate_generator = PrimitiveInternalCoordinates(geo_molecule)

    internal_coordinate_types = {
        "Bond": Distance,
        "Angle": Angle,
        "Dihedral": Dihedral,
        "Improper": OutOfPlane,
    }

    internal_coordinates = {
        label: [
            (
                internal_coordinate.value(qm_conformer),
                internal_coordinate.value(mm_conformer),
            )
            for internal_coordinate in internal_coordinate_generator.Internals
            if isinstance(internal_coordinate, internal_coordinate_class)
        ]
        for label, internal_coordinate_class in internal_coordinate_types.items()
    }

    internal_coordinate_rmsd = {}

    for ic_type, ic_values in internal_coordinates.items():

        if len(ic_values) == 0:
            continue

        qm_values, mm_values = zip(*ic_values)

        qm_values = numpy.array(qm_values)
        mm_values = numpy.array(mm_values)

        # Converting from radians to degrees
        if ic_type in ["Angle", "Dihedral", "Improper"]:
            rmsd = compute_rmsd(qm_values*180/numpy.pi, mm_values*180/numpy.pi, 360)
        else:
            rmsd = compute_rmsd(qm_values, mm_values)

        internal_coordinate_rmsd[ic_type] = float(rmsd)

    fb_objective = 0.0
    for key, value in internal_coordinate_rmsd.items():
        if key == "Bond":
            fb_objective += (1/0.1)*internal_coordinate_rmsd["Bond"]
        elif key == "Angle":
            fb_objective += (1/8.0)*internal_coordinate_rmsd["Angle"]
        elif key == "Dihedral":
            fb_objective += (1/20.0) * internal_coordinate_rmsd["Dihedral"]
        elif key == "Improper":
            fb_objective += (1/20.0) * internal_coordinate_rmsd["Improper"]

    return internal_coordinate_rmsd, fb_objective


def _compute_metrics(
    label: str, input_dictionary: Dict[str, Any], index: int
) -> pandas.DataFrame:

    record_id_tag = input_dictionary["record-id-tag"]

    qm_tag = input_dictionary["qm-tag"]
    mm_tag = input_dictionary["mm-tag"]

    qm_molecule_stream = oechem.oemolistream(
        input_dictionary["qm-structure"] + f"-{index}.sdf"
    )
    qm_molecule_stream.SetConfTest(oechem.OEAbsoluteConfTest(False))

    mm_molecule_stream = oechem.oemolistream(
        input_dictionary["mm-structures"][label] + f"-{index}.sdf"
    )
    mm_molecule_stream.SetConfTest(oechem.OEAbsoluteConfTest(False))

    metrics = []

    for oe_qm_molecule, oe_mm_molecule in zip(
        qm_molecule_stream.GetOEMols(), mm_molecule_stream.GetOEMols()
    ):

        conformer_metrics = []

        for conformer_index, (oe_qm_conformer, oe_mm_conformer) in enumerate(
            zip(oe_qm_molecule.GetConfs(), oe_mm_molecule.GetConfs())
        ):

            smiles = oechem.OEGetSDData(oe_qm_conformer, "SMILES QCArchive")
            record_id = oechem.OEGetSDData(oe_qm_conformer, record_id_tag)

            assert (
                oechem.OEGetSDData(oe_mm_conformer, record_id_tag) == record_id
            ), f"conformer mismatch {oe_qm_conformer.GetTitle()}"

            qm_energy = float(oechem.OEGetSDData(oe_qm_conformer, qm_tag))
            mm_energy = float(oechem.OEGetSDData(oe_mm_conformer, mm_tag))

            rmsd = oechem.OERMSD(
                oe_qm_conformer,
                oe_mm_conformer,
                RMSD_AUTOMORPH,
                RMSD_HEAVY_ONLY,
                RMSD_OVERLAY,
            )

            qm_molecule = Molecule.from_openeye(
                oechem.OEMol(oe_qm_conformer), allow_undefined_stereo=True
            )
            mm_molecule = Molecule.from_openeye(
                oechem.OEMol(oe_mm_conformer), allow_undefined_stereo=True
            )

            internal_coordinate_rmsds, fb_objective = _compute_internal_coordinate_rmsd(
                qm_molecule, qm_molecule.conformers[0], mm_molecule.conformers[0]
            )

            rd_qm_molecule = qm_molecule.to_rdkit()
            rd_mm_molecule = mm_molecule.to_rdkit()

            try:
                tfd = TorsionFingerprints.GetTFDBetweenMolecules(
                    rd_qm_molecule, rd_mm_molecule
                )

            except IndexError:

                print(
                    f"error calculating TFD for id={record_id} - possibly no "
                    f"non-terminal rotatable bonds found."
                )

                tfd = np.nan

            conformer_metrics.append(
                {
                    "SMILES": smiles,
                    "Conformer Idx": conformer_index,
                    "QM Energy": qm_energy,
                    "MM Energy": mm_energy,
                    "RMSD": rmsd,
                    **{
                        f"{ic_type} RMSD": ic_rmsd
                        for ic_type, ic_rmsd in internal_coordinate_rmsds.items()
                    },
                    "TDF": tfd,
                    "FB OBJECTIVE" : fb_objective,
                    "Record ID": record_id,
                }
            )

        # Compute the dde metric.
        qm_energies = numpy.array([metric["QM Energy"] for metric in conformer_metrics])
        mm_energies = numpy.array([metric["MM Energy"] for metric in conformer_metrics])

        lowest_qm_energy_idx = qm_energies.argmin()

        relative_qm_energies = qm_energies - qm_energies[lowest_qm_energy_idx]
        relative_mm_energies = mm_energies - mm_energies[lowest_qm_energy_idx]

        relative_qm_energies[lowest_qm_energy_idx] = numpy.nan
        relative_mm_energies[lowest_qm_energy_idx] = numpy.nan

        dde = relative_mm_energies - relative_qm_energies

        metrics.extend(
            {
                "Force Field": label,
                **{
                    key: value
                    for key, value in conformer_metric.items()
                    if "Energy" not in key
                },
                "ddE": conformer_dde,
            }
            for conformer_dde, conformer_metric in zip(dde, conformer_metrics)
        )

    return pandas.DataFrame(metrics)


@click.command()
@click.option(
    "--input",
    "input_path",
    type=click.Path(exists=True, dir_okay=False),
    default="03-force-fields.json",
)
@click.option(
    "--index",
    "index",
    type=click.INT,
)
@click.option(
    "--output",
    "output_path",
    type=click.Path(exists=False, dir_okay=False),
)
def main(input_path, index, output_path):

    with open(input_path) as file:
        input_dictionary = json.load(file)

    method_labels = [*input_dictionary["mm-structures"]]

    print(method_labels)
    metrics = pandas.concat(
        [_compute_metrics(label, input_dictionary, index) for label in method_labels]
    )
    metrics.to_csv(output_path, index=False)


if __name__ == "__main__":
    main()
