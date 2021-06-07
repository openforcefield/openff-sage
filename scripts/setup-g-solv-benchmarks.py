import requests_mock
from nonbonded.library.config import settings
from nonbonded.library.factories.inputs.benchmark import BenchmarkInputFactory
from nonbonded.library.models.authors import Author
from nonbonded.library.models.datasets import DataSet
from nonbonded.library.models.forcefield import ForceField
from nonbonded.library.models.projects import Benchmark
from nonbonded.library.models.results import OptimizationResult
from nonbonded.library.utilities.environments import ChemicalEnvironment
from openff.evaluator.forcefield import TLeapForceFieldSource
from openff.toolkit.typing.engines.smirnoff.forcefield import (
    ForceField as OFFForceField,
)

PROJECT_ID = "openff-force-fields"
STUDY_ID = "sage"

AUTHORS = [
    Author(
        name="Simon Boothroyd",
        email="simon.boothroyd@hotmail.com",
        institute="Boothroyd Scientific Consulting",
    ),
]

CHEMICAL_ENVIRONMENTS = [
    ChemicalEnvironment.Aqueous,
    ChemicalEnvironment.Alkane,
    ChemicalEnvironment.Alkene,
    ChemicalEnvironment.Alcohol,
    ChemicalEnvironment.CarbonylHydrate,
    ChemicalEnvironment.Hemiacetal,
    ChemicalEnvironment.Acetal,
    ChemicalEnvironment.Hemiaminal,
    ChemicalEnvironment.Aminal,
    ChemicalEnvironment.Thioacetal,
    ChemicalEnvironment.CarboxylicAcidEster,
    ChemicalEnvironment.Ether,
    ChemicalEnvironment.Aldehyde,
    ChemicalEnvironment.Ketone,
    ChemicalEnvironment.Aromatic,
    ChemicalEnvironment.CarboxylicAcidPrimaryAmide,
    ChemicalEnvironment.CarboxylicAcidSecondaryAmide,
    ChemicalEnvironment.CarboxylicAcidTertiaryAmide,
    ChemicalEnvironment.PrimaryAmine,
    ChemicalEnvironment.SecondaryAmine,
    ChemicalEnvironment.TertiaryAmine,
    ChemicalEnvironment.Cyanate,
    ChemicalEnvironment.Isocyanate,
    ChemicalEnvironment.Heterocycle,
    ChemicalEnvironment.AlkylFluoride,
    ChemicalEnvironment.ArylFluoride,
    ChemicalEnvironment.AlkylChloride,
    ChemicalEnvironment.ArylChloride,
    ChemicalEnvironment.AlkylBromide,
    ChemicalEnvironment.ArylBromide,
    ChemicalEnvironment.Thiol,
    ChemicalEnvironment.Thioaldehyde,
    ChemicalEnvironment.Thioketone,
    ChemicalEnvironment.Thioether,
    ChemicalEnvironment.Disulfide,
    ChemicalEnvironment.Thiourea,
    ChemicalEnvironment.Thiocyanate,
    ChemicalEnvironment.Isothiocyanate,
    ChemicalEnvironment.ThiocarboxylicAcid,
    ChemicalEnvironment.ThiocarboxylicAcidEster,
    ChemicalEnvironment.SulfonicAcidEster,
    ChemicalEnvironment.Sulfone,
    ChemicalEnvironment.PhosphonicAcid,
    ChemicalEnvironment.PhosphoricAcid,
    ChemicalEnvironment.PhosphoricAcidEster,
]


def setup_benchmark(
    benchmark_id,
    name,
    description,
    optimization_id,
    force_field,
    data_set_id,
):

    if force_field is not None:
        force_field = ForceField.from_openff(force_field)

    benchmark = Benchmark(
        project_id=PROJECT_ID,
        study_id=STUDY_ID,
        id=benchmark_id,
        name=name,
        description=description,
        optimization_id=optimization_id,
        force_field=force_field,
        analysis_environments=CHEMICAL_ENVIRONMENTS,
        test_set_ids=[data_set_id],
    )

    data_set_path = (
        f"../data-set-curation/physical-property/benchmarks/"
        f"data-sets/{data_set_id}.json"
    )

    test_set = DataSet.parse_file(data_set_path)

    for i, entry in enumerate(test_set.entries):
        entry.id = i + 1

    results = None

    if optimization_id is not None:

        results = OptimizationResult.from_rest(
            project_id=PROJECT_ID, study_id=STUDY_ID, model_id=benchmark.optimization_id
        )

    with requests_mock.Mocker() as request_mock:

        request_mock.get(
            f"{settings.API_URL}/datasets/phys-prop/{test_set.id}", text=test_set.json()
        )

        if results is not None:

            request_mock.get(
                (
                    f"https://nonbonded.herokuapp.com/api/dev/"
                    f"projects/{results.project_id}/"
                    f"studies/{results.study_id}/"
                    f"optimizations/{results.id}/results/"
                ),
                text=results.json(),
            )

        BenchmarkInputFactory.generate(
            benchmark,
            conda_environment="openff-force-fields",
            max_time="168:00",
            evaluator_preset="lilac-dask",
            evaluator_port=8000,
            n_evaluator_workers=70,
            include_results=False,
        )


def main():

    setup_benchmark(
        benchmark_id="mnsol-gaff-2-11",
        name="MNSol Benchmark of GAFF 2.11 + AM1BCC",
        description=(
            "An benchmark of the GAFF 2.11 force field + the AM1BCC charge model "
            "against a subset of the MNSol data set."
        ),
        optimization_id=None,
        force_field=TLeapForceFieldSource("leaprc.gaff2"),
        data_set_id="sage-mnsol-test-v1",
    )
    setup_benchmark(
        benchmark_id="mnsol-openff-1-3-0",
        name="MNSol Benchmark of OpenFF 1.3.0",
        description=(
            "An benchmark of the OpenFF 1.3.0 force field against a subset of the "
            "MNSol data set."
        ),
        optimization_id=None,
        force_field=OFFForceField("openff-1.3.0.offxml"),
        data_set_id="sage-mnsol-test-v1"
    )
    setup_benchmark(
        benchmark_id="mnsol-vdw-v1",
        name="MNSol Benchmark of VdW Parameters V1",
        description=(
            "An benchmark of the force field produced by the `vdw-v1` "
            "optimization against a subset of the MNSol data set."
        ),
        optimization_id="vdw-v1",
        force_field=None,
        data_set_id="sage-mnsol-test-v1"
    )

    setup_benchmark(
        benchmark_id="fsolv-gaff-2-11",
        name="FreeSolv Benchmark of GAFF 2.11 + AM1BCC",
        description=(
            "An benchmark of the GAFF 2.11 force field + the AM1BCC charge model "
            "against a subset of the FreeSolv data set."
        ),
        optimization_id=None,
        force_field=TLeapForceFieldSource("leaprc.gaff2"),
        data_set_id="sage-fsolv-test-v1",
    )
    setup_benchmark(
        benchmark_id="fsolv-openff-1-3-0",
        name="FreeSolv Benchmark of OpenFF 1.3.0",
        description=(
            "An benchmark of the OpenFF 1.3.0 force field against a subset of the "
            "FreeSolv data set."
        ),
        optimization_id=None,
        force_field=OFFForceField("openff-1.3.0.offxml"),
        data_set_id="sage-fsolv-test-v1",
    )
    setup_benchmark(
        benchmark_id="fsolv-vdw-v1",
        name="FreeSolv Benchmark of VdW Parameters V1",
        description=(
            "An benchmark of the force field produced by the `vdw-v1` "
            "optimization against a subset of the FreeSolv data set."
        ),
        optimization_id="vdw-v1",
        force_field=None,
        data_set_id="sage-fsolv-test-v1",
    )


if __name__ == "__main__":
    main()
