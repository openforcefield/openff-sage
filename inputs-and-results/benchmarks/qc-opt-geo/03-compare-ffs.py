import functools
import json
import pickle
from multiprocessing import Pool
from typing import Any, Dict, Tuple

import click
import numpy
import numpy as np
import openeye.oechem as oechem
from openff.toolkit.topology import Molecule
from rdkit.Chem import TorsionFingerprints

# Define the RMSD calculation parameters
RMSD_AUTOMORPH = True  # take into acct symmetry related transformations
RMSD_HEAVY_ONLY = False  # do consider hydrogen atoms for automorphisms
RMSD_OVERLAY = True  # find the lowest possible RMSD


def _compare_force_field(
    label: str, input_dictionary: Dict[str, Any]
) -> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]:

    record_id_tag = input_dictionary["record-id-tag"]

    qm_tag = input_dictionary["qm-tag"]
    mm_tag = input_dictionary["mm-tag"]

    qm_molecule_stream = oechem.oemolistream(input_dictionary["qm-structure"])
    qm_molecule_stream.SetConfTest(oechem.OEAbsoluteConfTest(False))
    mm_molecule_stream = oechem.oemolistream(input_dictionary["mm-structures"][label])
    mm_molecule_stream.SetConfTest(oechem.OEAbsoluteConfTest(False))

    full_energies = []
    full_rmsds = []
    full_tfds = []
    full_smiles = []

    for qm_molecule, mm_molecule in zip(
        qm_molecule_stream.GetOEMols(), mm_molecule_stream.GetOEMols()
    ):

        qm_energies = []
        mm_energies = []

        conformer_rmsds = []
        conformer_tfds = []

        conformer_smiles = []

        for qm_conformer, mm_conformer in zip(
            qm_molecule.GetConfs(), mm_molecule.GetConfs()
        ):

            conformer_smiles.append(
                oechem.OEGetSDData(qm_conformer, "SMILES QCArchive")
            )

            record_id = oechem.OEGetSDData(qm_conformer, record_id_tag)

            assert (
                oechem.OEGetSDData(mm_conformer, record_id_tag) == record_id
            ), f"conformer mismatch {qm_conformer.GetTitle()}"

            qm_energies.append(float(oechem.OEGetSDData(qm_conformer, qm_tag)))
            mm_energies.append(float(oechem.OEGetSDData(mm_conformer, mm_tag)))

            conformer_rmsds.append(
                oechem.OERMSD(
                    qm_conformer,
                    mm_conformer,
                    RMSD_AUTOMORPH,
                    RMSD_HEAVY_ONLY,
                    RMSD_OVERLAY,
                )
            )

            rd_qm_molecule = Molecule.from_openeye(
                oechem.OEMol(qm_conformer), allow_undefined_stereo=True
            ).to_rdkit()
            rd_mm_molecule = Molecule.from_openeye(
                oechem.OEMol(mm_conformer), allow_undefined_stereo=True
            ).to_rdkit()

            print(rd_qm_molecule.GetNumConformers())
            print(rd_mm_molecule.GetNumConformers())

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

            conformer_tfds.append(tfd)

        lowest_qm_energy_idx = qm_energies.index(min(qm_energies))

        relative_qm_energies = np.array(qm_energies) - qm_energies[lowest_qm_energy_idx]
        relative_mm_energies = np.array(mm_energies) - mm_energies[lowest_qm_energy_idx]

        relative_qm_energies = np.delete(relative_qm_energies, lowest_qm_energy_idx)
        relative_mm_energies = np.delete(relative_mm_energies, lowest_qm_energy_idx)

        conformer_rmsds.pop(lowest_qm_energy_idx)
        conformer_tfds.pop(lowest_qm_energy_idx)
        conformer_smiles.pop(lowest_qm_energy_idx)

        # subtract them to get ddE = dE (query method) - dE (ref method)
        conformer_energies = np.array(relative_mm_energies) - np.array(
            relative_qm_energies
        )

        full_energies.extend(conformer_energies.tolist())
        full_rmsds.extend(conformer_rmsds)
        full_tfds.extend(conformer_tfds)
        full_smiles.extend(conformer_smiles)

    return (
        numpy.array(full_energies),
        numpy.array(full_rmsds),
        numpy.array(full_tfds),
        numpy.array(full_smiles),
    )


@click.command()
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(exists=True, dir_okay=False),
    default="03-force-fields.json",
)
@click.option(
    "-o",
    "--output",
    "output_path",
    type=click.Path(exists=False, dir_okay=False),
    default="03-metrics.pkl",
)
def main(input_path, output_path):

    with open(input_path) as file:
        input_dictionary = json.load(file)

    method_labels = [*input_dictionary["mm-structures"]]

    with Pool(processes=len(method_labels)) as pool:

        metrics = list(
            pool.map(
                functools.partial(
                    _compare_force_field, input_dictionary=input_dictionary
                ),
                method_labels,
            )
        )

    metrics_per_label = {label: metric for label, metric in zip(method_labels, metrics)}

    with open(output_path, "wb") as file:
        pickle.dump(metrics_per_label, file)


if __name__ == "__main__":
    main()
