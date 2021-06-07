import os
from collections import defaultdict
from typing import List, Tuple

import click
import numpy
import pandas
import seaborn
from matplotlib import pyplot
from nonbonded.library.utilities.checkmol import analyse_functional_groups
from nonbonded.library.utilities.environments import ChemicalEnvironment
from tqdm import tqdm

seaborn.set_palette("colorblind")

ANALYSIS_ENVIRONMENTS = [
    ChemicalEnvironment.Aldehyde,
    ChemicalEnvironment.Ketene,
    ChemicalEnvironment.Alcohol,
    ChemicalEnvironment.Ether,
    ChemicalEnvironment.Disulfide,
    ChemicalEnvironment.PrimaryAmine,
    ChemicalEnvironment.SecondaryAmine,
    ChemicalEnvironment.TertiaryAmine,
    ChemicalEnvironment.AlkylChloride,
    ChemicalEnvironment.ArylChloride,
    ChemicalEnvironment.AlkylFluoride,
    ChemicalEnvironment.ArylFluoride,
    ChemicalEnvironment.CarboxylicAcid,
    ChemicalEnvironment.CarbamicAcidEster,
    ChemicalEnvironment.CarboxylicAcidPrimaryAmide,
    ChemicalEnvironment.CarboxylicAcidSecondaryAmide,
    ChemicalEnvironment.CarboxylicAcidTertiaryAmide,
    ChemicalEnvironment.Nitrile,
    ChemicalEnvironment.CarbamicAcid,
    ChemicalEnvironment.Urea,
    ChemicalEnvironment.NitroCompound,
    ChemicalEnvironment.Sulfonamide,
    ChemicalEnvironment.SulfuricAcidDeriv,
]


def draw_step_plot(
    plot_data: List[numpy.ndarray],
    method_labels: List[str],
    metric: str,
    x_range: Tuple[float, float],
    output_path: str,
    n_bins: int = 16,
    bootstrap_iterations: int = 2000,
    percentile=0.95,
):
    """Draw a step plot of data segregated by each method.

    Parameters
    ----------
    plot_data
        A list of the arrays where ``plot_data[i]`` contains the particular metric of
        interest associated with ``method_labels[i]``.
    method_labels
        The label associated with each data series.
    metric
        The metric contained in the ``plot_data``.
    x_range
        A tuple of the min and max values to use for the x axis.
    output_path
        The path to save the final plot to.
    n_bins
        The number of bins to split the data into.
    bootstrap_iterations
        The number of bootstrap iterations to perform when computing the error bars.
    percentile
        The percentile of the confidence interval to display.
    """

    sample_count = len(plot_data[0])

    sample_densities = {
        method: numpy.zeros((bootstrap_iterations, n_bins)) for method in method_labels
    }

    if x_range is None:
        x_range = (
            min(numpy.nanmin(data) for data in plot_data),
            max(numpy.nanmax(data) for data in plot_data),
        )

    for sample_index in range(bootstrap_iterations):

        samples_indices = numpy.random.randint(
            low=0, high=sample_count, size=sample_count
        )

        for method, method_data in zip(method_labels, plot_data):

            sample_data = method_data[samples_indices]
            sample_density, _ = numpy.histogram(
                sample_data, bins=n_bins, range=x_range, density=False
            )

            sample_densities[method][sample_index] = sample_density

    lower_percentile_index = int(bootstrap_iterations * (1 - percentile) / 2)
    upper_percentile_index = int(bootstrap_iterations * (1 + percentile) / 2)

    confidence_intervals = {}

    error_bar_x_values = {}
    error_bar_y_values = {}

    for method, method_data in zip(method_labels, plot_data):

        sorted_samples = numpy.sort(sample_densities[method], axis=0)

        confidence_intervals[method] = (
            sorted_samples[lower_percentile_index],
            sorted_samples[upper_percentile_index],
        )

        error_bar_x_values[method] = (
            numpy.arange(n_bins) / n_bins * (x_range[1] - x_range[0])
        ) + x_range[0]
        error_bar_y_values[method], _ = numpy.histogram(
            method_data, bins=n_bins, range=x_range, density=False
        )

    plot_frame = pandas.DataFrame(
        [
            {"method": method_label, metric: data_point}
            for data_series, method_label in zip(plot_data, method_labels)
            for data_point in data_series
        ]
    )

    plot = seaborn.displot(
        data=plot_frame,
        x=metric,
        hue="method",
        kind="hist",
        stat="count",
        bins=n_bins,
        binrange=x_range,
        element="step",
        fill=False,
        aspect=2.5,
        height=3,
        lw=1.25,
        facet_kws={"legend_out": False},
    )

    for i, method in enumerate(method_labels):

        x_offset = (x_range[1] - x_range[0]) / n_bins / len(method_labels) * (i + 0.5)

        pyplot.errorbar(
            error_bar_x_values[method] + x_offset,
            error_bar_y_values[method],
            yerr=numpy.array(
                [
                    numpy.abs(
                        confidence_intervals[method][0] - error_bar_y_values[method]
                    ),
                    numpy.abs(
                        confidence_intervals[method][1] - error_bar_y_values[method]
                    ),
                ]
            ),
            marker=None,
            linestyle="none",
        )

    # Remove the legend title.
    plot.legend.set_title(None)

    # Draw a vertical line at x=0 for visual reference
    pyplot.axvline(x=0.0, lw=1.5, ls="--", color="gray", clip_on=False)

    pyplot.gcf().set_size_inches(7, 3)

    # adjust font sizes
    pyplot.xlabel(metric, fontsize=14)
    pyplot.ylabel("Density", fontsize=14)

    # save with transparency for overlapping plots
    pyplot.savefig(output_path, transparent=True, bbox_inches="tight")
    pyplot.close("all")


def draw_box_plot(
    plot_data: List[numpy.ndarray],
    method_labels: List[str],
    metric: str,
    y_range: Tuple[float, float],
    output_path: str,
):
    """Draw a violin plot of data segregated by each method.

    Parameters
    ----------
    plot_data
        A list of the arrays where ``plot_data[i]`` contains the particular metric of
        interest associated with ``method_labels[i]``.
    method_labels
        The label associated with each data series.
    metric
        The metric contained in the ``plot_data``.
    output_path
        The path to save the final plot to.
    """

    plot_frame = pandas.DataFrame(
        [
            {"method": method_label, metric: data_point}
            for data_series, method_label in zip(plot_data, method_labels)
            for data_point in data_series
        ]
    )

    seaborn.catplot(data=plot_frame, x="method", y=metric, kind="box", dodge=False)

    pyplot.xticks(rotation=45)

    if y_range is not None:
        pyplot.ylim(y_range)

    pyplot.xlabel("Force Field", fontsize=14)
    pyplot.ylabel(metric, fontsize=14)

    # save with transparency for overlapping plots
    pyplot.savefig(output_path, transparent=True, bbox_inches="tight")
    pyplot.close("all")


def draw_plots(energies, rmsds, tfds, method_labels, output_directory):
    """Draw step and violin plots of a set of ddE, RMSD, and TFD metrics."""

    os.makedirs(output_directory, exist_ok=True)

    draw_box_plot(
        energies,
        method_labels,
        metric="ddE (kcal/mol)",
        y_range=(-15.0, 15.0),
        output_path=os.path.join(output_directory, "box-dde.png"),
    )
    draw_step_plot(
        energies,
        method_labels,
        metric="ddE (kcal/mol)",
        x_range=(-15.0, 15.0),
        output_path=os.path.join(output_directory, "step-dde.png"),
    )

    for key, values in rmsds.items():
        draw_box_plot(
            values,
            method_labels,
            metric=fr"{key} ($\mathrm{{\AA}}$)",
            y_range=(0, 3) if key == "RMSD" else None,
            output_path=os.path.join(
                output_directory, f"box-{key.lower().replace(' ', '_')}.png"
            ),
        )
        draw_step_plot(
            values,
            method_labels,
            metric=fr"{key} ($\mathrm{{\AA}}$)",
            x_range=(0, 3) if key == "RMSD" else None,
            output_path=os.path.join(
                output_directory, f"step-{key.lower().replace(' ', '_')}.png"
            ),
        )

    draw_box_plot(
        tfds,
        method_labels,
        metric="TFD",
        y_range=(0, 0.5),
        output_path=os.path.join(output_directory, "box-tfd.png"),
    )
    draw_step_plot(
        tfds,
        method_labels,
        metric="TFD",
        x_range=(0, 0.5),
        output_path=os.path.join(output_directory, "step-tfd.png"),
    )


@click.command()
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(exists=True, dir_okay=False),
    default="03-metrics.csv",
)
def main(input_path):

    print("1) Plotting full metrics")

    metrics_frame = pandas.read_csv(input_path)

    method_labels = sorted(metrics_frame["Force Field"].unique())

    metrics = []

    for method_label in method_labels:

        method_frame = metrics_frame[metrics_frame["Force Field"] == method_label]

        metrics.append(
            (
                method_frame["SMILES"].values,
                method_frame["ddE"].values,
                method_frame["TDF"].values,
                method_frame["RMSD"].values,
                method_frame["Bond RMSD"].values,
                method_frame["Angle RMSD"].values,
                method_frame["Dihedral RMSD"].values,
                method_frame["Improper RMSD"].values,
            )
        )

    # Plot the full statistics
    (
        smiles,
        energies,
        tfds,
        full_rmsds,
        bond_rmsds,
        angle_rmsds,
        proper_rmsds,
        improper_rmsds,
    ) = zip(*metrics)

    draw_plots(
        energies,
        {
            "RMSD": full_rmsds,
            "Bond RMSD": bond_rmsds,
            "Angle RMSD": angle_rmsds,
            "Proper RMSD": proper_rmsds,
            "Improper RMSD": improper_rmsds,
        },
        tfds,
        method_labels,
        os.path.join("04-outputs", "full"),
    )

    # Split the data into per environment results
    per_environment_smiles = defaultdict(set)

    print("2) Splitting test molecules by chemical environment")

    for pattern in tqdm(smiles[0]):

        chemical_environments = analyse_functional_groups(pattern)

        if chemical_environments is None:
            print(f"skipping {pattern} - checkmol could not assign groups")

        for chemical_environment in chemical_environments:

            if chemical_environment not in ANALYSIS_ENVIRONMENTS:
                continue

            per_environment_smiles[chemical_environment].add(pattern)

    print("3) Splitting test data by chemical environment")

    per_environment_metrics = defaultdict(list)

    for method_smiles, method_energies, method_tfds, method_rmsds, *_ in metrics:

        method_data = pandas.DataFrame(
            {
                "Energy": method_energies,
                "RMSD": method_rmsds,
                "TFD": method_tfds,
                "SMILES": method_smiles,
            }
        )

        for chemical_environment, environment_smiles in per_environment_smiles.items():

            chemical_environment_data = method_data[
                method_data["SMILES"].isin(environment_smiles)
            ]

            per_environment_metrics[chemical_environment].append(
                (
                    chemical_environment_data["Energy"].values,
                    chemical_environment_data["RMSD"].values,
                    chemical_environment_data["TFD"].values,
                    chemical_environment_data["SMILES"].values,
                )
            )

    print("3) Plotting per environment metrics")

    for chemical_environment, environment_metrics in per_environment_metrics.items():

        environment_energies, environment_rmsds, environment_tfds, _ = zip(
            *environment_metrics
        )

        draw_plots(
            environment_energies,
            {"RMSD": environment_rmsds},
            environment_tfds,
            method_labels,
            os.path.join("04-outputs", chemical_environment.value.lower()),
        )


if __name__ == "__main__":
    main()
