import base64
import io
import json
import os
import urllib
from collections import defaultdict
from typing import Dict, Tuple

import click
import numpy
import pandas
import seaborn
from matplotlib import pyplot
from tqdm import tqdm

seaborn.set_palette("colorblind")


def _extract_molecule_profiles(
    energy_frame: pandas.DataFrame, smiles: str
) -> Tuple[numpy.ndarray, numpy.ndarray, Dict[str, numpy.ndarray]]:
    """Extracts the QM and MM energies profiles from a master pandas data frame
    for a particular molecule defined by its SMILES pattern.

    Returns
        The grid angles, QM profile and MM profiles for each force field.
    """

    smiles_frame = energy_frame[energy_frame["SMILES"] == smiles]

    qm_energies = {}
    mm_energies = defaultdict(dict)

    for _, row in smiles_frame.iterrows():
        qm_energies[int(row["Grid ID"])] = row["QM Energy"]
        mm_energies[row["Force Field"]][int(row["Grid ID"])] = row["MM Energy"]

    grid_ids = sorted(grid_id for grid_id in qm_energies)

    qm_energies = numpy.array([qm_energies[grid_id] for grid_id in grid_ids])

    mm_energies = {
        mm_label: numpy.array([mm_energies[mm_label][grid_id] for grid_id in grid_ids])
        for mm_label in mm_energies
    }

    return numpy.array(grid_ids), qm_energies, mm_energies


def _compute_profile_rmse(
    qm_energies: numpy.ndarray, mm_energies: numpy.ndarray
) -> float:
    """Compute the RMSE between a QM and an MM torsion profile after optimally
    superimposing the two."""

    shift = (qm_energies - mm_energies).mean()
    mm_energies = mm_energies + shift

    delta = numpy.sqrt(
        (qm_energies - mm_energies).dot(qm_energies - mm_energies) / len(qm_energies)
    )

    return float(delta)


def _smiles_to_img(smiles: str) -> str:
    """A helper function for creating an HTML image from a SMILES pattern."""

    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    rdkit_molecule = Chem.MolFromSmiles(smiles)

    tagged_bonds = [
        bond.GetIdx()
        for bond in rdkit_molecule.GetBonds()
        if bond.GetBeginAtom().GetAtomMapNum() != 0
        and bond.GetEndAtom().GetAtomMapNum() != 0
    ]

    # Generate a set of 2D coordinates.
    Chem.rdDepictor.Compute2DCoords(rdkit_molecule)

    drawer = rdMolDraw2D.MolDraw2DSVG(400, 300)
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        rdkit_molecule,
        highlightBonds=tagged_bonds,
    )
    drawer.FinishDrawing()

    svg_content = drawer.GetDrawingText()

    return (
        f'<img src="data:image/svg+xml;base64,'
        f'{base64.b64encode(svg_content.encode()).decode()}"/>'
    )


def plot_torsion_profile(
    angles: numpy.ndarray,
    qm_energies: numpy.ndarray,
    mm_energies: Dict[str, numpy.ndarray]
) -> str:
    """Plots a set of QM and MM torsion profiles and returns the result as a HTML
    encoded image.
    """

    figure = pyplot.figure(figsize=(6.0, 4.0))

    pyplot.plot(angles, qm_energies, label="QM")

    for mm_label in mm_energies:

        shift = (qm_energies - mm_energies[mm_label]).mean()
        mm_energies_method = mm_energies[mm_label] + shift

        pyplot.plot(angles, mm_energies_method, label=mm_label)

    pyplot.xlabel("Angle")
    pyplot.ylabel("Energy (kcal / mol)")

    pyplot.legend()

    pyplot.tight_layout()

    with io.BytesIO() as buffer:

        pyplot.gcf().savefig(buffer, format="png")
        buffer.seek(0)
        plot_string = base64.b64encode(buffer.read())

    pyplot.close(figure)

    img_html = f'<img src="data:image/png;base64,{urllib.parse.quote(plot_string)}"/>'
    return img_html


def plot_per_force_field_rmse(
    energy_frame: pandas.DataFrame,
    output_path: str,
    bootstrap_iterations: int = 2000,
    percentile=0.95,
):
    """Plots the mean RMSE between the QM torsion profile and the MM torsion profile
    computed for a set of force fields, showing bootstrapped error bars."""

    rmse_values = defaultdict(list)

    for smiles in tqdm(energy_frame["SMILES"].unique()):

        angles, qm_energies, mm_energies = _extract_molecule_profiles(
            energy_frame, smiles
        )

        for mm_label in mm_energies:

            rmse_values[mm_label].append(
                _compute_profile_rmse(
                    numpy.array(qm_energies), numpy.array(mm_energies[mm_label])
                )
            )

    method_labels = [*rmse_values]
    plot_data = [numpy.array(rmse_values[mm_label]) for mm_label in method_labels]

    # Bootstrap the data to compute the error bars.
    sample_count = len(plot_data[0])

    sample_values = {
        method: numpy.zeros((bootstrap_iterations, 1)) for method in method_labels
    }

    for sample_index in range(bootstrap_iterations):

        samples_indices = numpy.random.randint(
            low=0, high=sample_count, size=sample_count
        )

        for method, method_data in zip(method_labels, plot_data):
            sample_values[method][sample_index] = method_data[samples_indices].mean()

    lower_percentile_index = int(bootstrap_iterations * (1 - percentile) / 2)
    upper_percentile_index = int(bootstrap_iterations * (1 + percentile) / 2)

    # Produce the figure.
    figure = pyplot.figure(figsize=(4.0, 4.0))

    for i, (method, method_data) in enumerate(zip(method_labels, plot_data)):

        sorted_samples = numpy.sort(sample_values[method], axis=0)

        mean_value = method_data.mean()

        confidence_intervals = numpy.array(
            [
                [float(numpy.abs(mean_value - sorted_samples[lower_percentile_index]))],
                [float(numpy.abs(mean_value - sorted_samples[upper_percentile_index]))],
            ]
        )

        pyplot.errorbar([i], [mean_value], yerr=confidence_intervals, marker="x")

    pyplot.xticks(
        ticks=[*range(len(method_labels))], labels=method_labels, rotation=45.0
    )

    pyplot.tight_layout()
    pyplot.savefig(output_path)
    pyplot.close(figure)


def plot_per_molecule_profile(energy_frame: pandas.DataFrame, output_path: str):

    plot_rows = []

    for smiles in tqdm(energy_frame["SMILES"].unique()):

        angles, qm_energies, mm_energies = _extract_molecule_profiles(
            energy_frame, smiles
        )

        plot_rows.append(
            {
                "Molecule": _smiles_to_img(smiles),
                "Profile": plot_torsion_profile(
                    angles=angles, qm_energies=qm_energies, mm_energies=mm_energies
                ),
                **{
                    f"{mm_label} RMSE": _compute_profile_rmse(
                        numpy.array(qm_energies), numpy.array(mm_energies[mm_label])
                    )
                    for mm_label in mm_energies
                },
            }
        )

    plot_frame = pandas.DataFrame(plot_rows)

    with open(output_path, "w") as file:
        file.write(plot_frame.to_html(escape=False))


@click.command()
@click.option(
    "-i",
    "--input",
    "input_path",
    default="01-residuals.json",
    show_default=True,
    type=click.STRING,
    help="The path to the JSON file containing the energies and residuals produced by "
    "step 01.",
)
@click.option(
    "-o",
    "--output",
    "output_directory",
    default="02-plots",
    show_default=True,
    type=click.STRING,
    help="The directory to save the plots into.",
)
def main(input_path, output_directory):

    os.makedirs(output_directory, exist_ok=True)

    # Store all of the data in a pandas frame for easier manipulation.
    with open(input_path) as file:
        results = json.load(file)

    record_smiles = results["record_smiles"]
    energies = results["energies"]

    energy_frame = pandas.DataFrame(
        [
            {
                "QC Record": record_id,
                "SMILES": record_smiles[record_id],
                "Grid ID": grid_id,
                "Force Field": force_field_label,
                "QM Energy": energies[record_id][grid_id][force_field_label][0],
                "MM Energy": energies[record_id][grid_id][force_field_label][1],
                "MM Residual": energies[record_id][grid_id][force_field_label][2],
            }
            for record_id in energies
            for grid_id in energies[record_id]
            for force_field_label in energies[record_id][grid_id]
        ]
    )

    # Plot the per molecule profiles and global RMSE values.
    print("1) Plotting RMSE values")

    plot_per_force_field_rmse(energy_frame, os.path.join(output_directory, "rmse.png"))

    print("2) Plotting per molecule profiles")

    plot_per_molecule_profile(
        energy_frame, os.path.join(output_directory, "torsion-profiles.html")
    )


if __name__ == "__main__":
    main()
