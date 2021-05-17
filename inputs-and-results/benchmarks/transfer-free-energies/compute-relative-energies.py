import os
from collections import defaultdict
from typing import Dict, Optional, Tuple

import numpy
import pandas
import seaborn
from matplotlib import pyplot
from nonbonded.library.models.datasets import DataSet, DataSetEntry
from nonbonded.library.models.results import (
    BenchmarkResult,
    DataSetResult,
    DataSetResultEntry,
)
from nonbonded.library.plotting.utilities import sort_categories_key
from nonbonded.library.statistics.statistics import StatisticType, compute_statistics
from openff.evaluator import unit
from openff.toolkit.topology import Molecule


def _entries_to_measurements(
    data_entry: DataSetEntry, result_entry: DataSetResultEntry
) -> Tuple[unit.Measurement, unit.Measurement]:

    reference_measurement = unit.Measurement(
        data_entry.value, data_entry.std_error, data_entry.units
    )
    estimated_measurement = unit.Measurement(
        result_entry.estimated_value, result_entry.estimated_std_error, data_entry.units
    )

    return reference_measurement, estimated_measurement


def _partition_results_by_solute(
    reference_data_set: DataSet, results: DataSetResult
) -> Dict[str, Dict[str, Tuple[DataSetEntry, DataSetResultEntry]]]:

    reference_entries = {entry.id: entry for entry in reference_data_set.entries}

    per_solute_results = defaultdict(dict)

    # Partition the results by solute. We re-create the SMILES patterns just incase
    # different toolkits were used to generate them.
    for result_entry in results.result_entries:

        reference_entry = reference_entries[result_entry.reference_id]

        solute = next(filter(lambda x: x.role == "Solute", reference_entry.components))
        solvent = next(
            filter(lambda x: x.role == "Solvent", reference_entry.components)
        )

        solute_smiles = Molecule.from_smiles(solute.smiles).to_smiles()
        solvent_smiles = Molecule.from_smiles(solvent.smiles).to_smiles()

        per_solute_results[solute_smiles][solvent_smiles] = (
            reference_entry,
            result_entry,
        )

    return per_solute_results


def _to_transfer_free_energy_entries(
    reference_entry_1: DataSetEntry,
    estimated_entry_1: DataSetResultEntry,
    reference_entry_2: DataSetEntry,
    estimated_entry_2: DataSetResultEntry,
) -> pandas.DataFrame:

    reference_g_solv_1, estimated_g_solv_1 = _entries_to_measurements(
        reference_entry_1, estimated_entry_1
    )
    reference_g_solv_2, estimated_g_solv_2 = _entries_to_measurements(
        reference_entry_2, estimated_entry_2
    )

    reference_transfer_free_energy = reference_g_solv_2 - reference_g_solv_1
    estimated_transfer_free_energy = estimated_g_solv_2 - estimated_g_solv_1

    categories = set()

    for category_1 in estimated_entry_1.categories:

        solute_category_1 = category_1.split(" + ")[1]
        solvent_category_1 = category_1.split(" + ")[0]

        for category_2 in estimated_entry_2.categories:

            solute_category_2 = category_2.split(" + ")[1]
            solvent_category_2 = category_2.split(" + ")[0]

            categories.add(
                f"{solute_category_1} + {solvent_category_1} -> {solvent_category_2}"
            )
            categories.add(
                f"{solute_category_2} + {solvent_category_1} -> {solvent_category_2}"
            )

    categories = [
        category.replace(" (x=1.0)", "").replace(" (n=1)", "").replace(" + ", " - ")
        for category in categories
    ]

    return pandas.DataFrame(
        [
            {
                "Reference Value": reference_transfer_free_energy.magnitude.n,
                "Reference Std": reference_transfer_free_energy.magnitude.s,
                "Estimated Value": estimated_transfer_free_energy.magnitude.n,
                "Estimated Std": estimated_transfer_free_energy.magnitude.s,
                "Category": category,
            }
            for category in categories
        ]
    )


def _compute_transfer_free_energies(
    per_solute_results, solvent_1: str, solvent_2: str
) -> pandas.DataFrame:

    result_frames = []

    for solute in per_solute_results:

        if solvent_1 not in per_solute_results[solute]:

            print(
                f"skipping {solute} + {solvent_1} -> "
                f"{solvent_2} - {solute} + {solvent_1} missing"
            )
            continue

        if solvent_2 not in per_solute_results[solute]:

            print(
                f"skipping {solute} + {solvent_1} -> "
                f"{solvent_2} - {solute} + {solvent_2} missing"
            )
            continue

        reference_entry_1, estimated_entry_1 = per_solute_results[solute][solvent_1]
        reference_entry_2, estimated_entry_2 = per_solute_results[solute][solvent_2]

        result_frames.append(
            _to_transfer_free_energy_entries(
                reference_entry_1,
                estimated_entry_1,
                reference_entry_2,
                estimated_entry_2,
            )
        )

    if len(result_frames) == 0:
        return pandas.DataFrame()

    return pandas.concat(result_frames, ignore_index=True, sort=False)


def _results_frame_to_rmsd_frame(
    results_frame: pandas.DataFrame,
    category: Optional[str],
    bootstrap_iterations: int,
) -> pandas.DataFrame:

    if category is not None:
        results_frame = results_frame[results_frame["Category"] == category]

    bulk_statistics, _, bulk_statistics_ci = compute_statistics(
        measured_values=results_frame["Reference Value"].values,
        measured_std=results_frame["Reference Std"].values,
        estimated_values=results_frame["Estimated Value"].values,
        estimated_std=results_frame["Estimated Std"].values,
        bootstrap_iterations=bootstrap_iterations,
        statistic_types=[StatisticType.RMSE],
    )

    statistics_rows = []

    for statistic_type in bulk_statistics:

        statistics_rows.append(
            {
                "Category": category,
                "Value": bulk_statistics[statistic_type],
                "Lower 95 CI": bulk_statistics_ci[statistic_type][0],
                "Upper 95 CI": bulk_statistics_ci[statistic_type][1],
            }
        )

    return pandas.DataFrame(statistics_rows)


def _plot_per_category_rmsd(rmsd_frame: pandas.DataFrame):

    force_field_names = sorted(rmsd_frame["Force Field"].unique())

    rmsd_frame = rmsd_frame[rmsd_frame["Category"].notna()]

    categories = sorted(
        rmsd_frame["Category"].unique(), key=sort_categories_key, reverse=True
    )
    category_indices = [x * 2 for x in range(1, len(categories) + 1)]

    n_categories = len(categories)

    with seaborn.axes_style("white"):
        figure, axis = pyplot.subplots(figsize=(10.0, 1.5 + (0.5 * n_categories)))

        shifts = numpy.linspace(-0.5, 0.5, len(force_field_names))

        # Plot the data for each force field.
        for index, force_field_name in enumerate(force_field_names):

            benchmark_plot_data = rmsd_frame[
                rmsd_frame["Force Field"] == force_field_name
            ]

            # Plot the error bars and with the circle marker on top.
            axis.barh(
                benchmark_plot_data["Category"]
                .replace(categories, category_indices)
                .values
                - shifts[index],
                benchmark_plot_data["Value"],
                xerr=[
                    numpy.abs(
                        benchmark_plot_data["Lower 95 CI"]
                        - benchmark_plot_data["Value"]
                    ),
                    numpy.abs(
                        benchmark_plot_data["Upper 95 CI"]
                        - benchmark_plot_data["Value"]
                    ),
                ],
                height=1.0 / max(len(force_field_names) - 1, 1),
                label=force_field_name,
            )

        # Add a simple legend.
        axis.legend()

        # Add title and axis names
        axis.set_yticks(category_indices)
        axis.set_yticklabels(categories)

        axis.set_xlim(left=0.0)
        axis.set_ylim(0.5, (len(categories) + 1) * 2.0 - 0.5)

        axis.set_xlabel("RMSD | kJ / mol")

    # Save the figure.
    figure.tight_layout()

    pyplot.savefig(os.path.join("plots-g-solv", "ddg-solv-na-aq-rmsd-cat.png"))
    pyplot.close(figure)


def _plot_bulk_rmsd(rsmd_frame: pandas.DataFrame):

    plot_data = rsmd_frame[rsmd_frame["Category"].isna()]
    force_field_names = sorted(plot_data["Force Field"].unique())

    plot = seaborn.FacetGrid(
        data=plot_data,
        hue="Force Field",
        sharey=False,
        height=4.0,
        aspect=0.6,
    )

    axis = plot.facet_axis(0, 0)

    # Plot the error bars and marker for each benchmark.
    for index, force_field_name in enumerate(force_field_names):

        benchmark_plot_data = plot_data[plot_data["Force Field"] == force_field_name]

        axis.errorbar(
            index,
            benchmark_plot_data["Value"],
            yerr=[
                benchmark_plot_data["Lower 95 CI"] - benchmark_plot_data["Value"],
                benchmark_plot_data["Upper 95 CI"] - benchmark_plot_data["Value"],
            ],
            label=force_field_name,
            capsize=5.0,
            capthick=1.5,
            linestyle="none",
            marker="o",
            markersize=10.0,
            markerfacecolor="none",
            color=seaborn.color_palette(n_colors=1)[0],
        )

    # Add axis names
    axis.set_xticks(range(0, len(force_field_names)))
    axis.set_xticklabels(force_field_names, rotation=90)
    axis.set_xlim(-0.5, len(force_field_names) - 0.5)

    axis.set_ylabel("RMSD | kJ / mol")

    plot.savefig(
        os.path.join("plots-g-solv", "ddg-solv-na-aq-rmsd.png"), bbox_inches="tight"
    )
    pyplot.close(plot.fig)


def main():

    force_field_names = ["openff-1-3-0", "vdw-v1"]

    data_set_directory = os.path.join(
        os.path.pardir,
        os.path.pardir,
        os.path.pardir,
        "data-set-curation",
        "physical-property",
        "benchmarks",
        "data-sets",
    )

    fsolv_test_set = DataSet.parse_file(
        os.path.join(data_set_directory, "sage-fsolv-test-v1.json")
    )
    mnsol_test_set = DataSet.parse_file(
        os.path.join(data_set_directory, "sage-mnsol-test-v1.json")
    )

    bulk_rmsd_frames = []

    for force_field_name in force_field_names:

        fsolv_results = BenchmarkResult.parse_file(
            os.path.join(
                f"fsolv-{force_field_name}", "analysis", "benchmark-results.json"
            )
        ).data_set_result

        mnsol_results = BenchmarkResult.parse_file(
            os.path.join(
                f"mnsol-{force_field_name}", "analysis", "benchmark-results.json"
            )
        ).data_set_result

        all_solvents = set()
        per_solute_results = defaultdict(dict)

        for test_set, results in [
            (fsolv_test_set, fsolv_results),
            (mnsol_test_set, mnsol_results),
        ]:

            for solute, solvents in _partition_results_by_solute(
                test_set, results
            ).items():

                all_solvents.update({*solvents})
                per_solute_results[solute].update(solvents)

        force_field_results_frame = pandas.concat(
            [
                _compute_transfer_free_energies(per_solute_results, solvent, "[H]O[H]")
                for solvent in all_solvents
                if solvent != "[H]O[H]"
            ],
            ignore_index=True,
            sort=False,
        )

        for category in [None] + [*force_field_results_frame["Category"].unique()]:

            force_field_rmsd_frame = _results_frame_to_rmsd_frame(
                force_field_results_frame, category, 1000
            )
            force_field_rmsd_frame["Force Field"] = force_field_name

            bulk_rmsd_frames.append(force_field_rmsd_frame)

    bulk_rmsd_frame = pandas.concat(bulk_rmsd_frames, ignore_index=True, sort=False)
    bulk_rmsd_frame.to_csv("ddg-solv-na-aq-rmsd.csv", index=False)

    bulk_rmsd_frame = pandas.read_csv("ddg-solv-na-aq-rmsd.csv")
    _plot_per_category_rmsd(bulk_rmsd_frame)
    _plot_bulk_rmsd(bulk_rmsd_frame)


if __name__ == "__main__":
    main()
