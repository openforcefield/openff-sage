import os
from typing import List

import pandas
from nonbonded.library.models.authors import Author
from nonbonded.library.models.datasets import DataSet
from openff.evaluator.datasets.curation.components import filtering, selection, thermoml
from openff.evaluator.datasets.curation.components.selection import State, TargetState
from openff.evaluator.datasets.curation.workflow import (
    CurationWorkflow,
    CurationWorkflowSchema,
)
from openff.evaluator.utils.checkmol import ChemicalEnvironment

AUTHORS = [
    Author(
        name="Simon Boothroyd",
        email="simon.boothroyd@hotmail.com",
        institute="Boothroyd Scientific Consulting Ltd.",
    ),
]

CHEMICAL_ENVIRONMENTS = [
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
    ChemicalEnvironment.Aqueous,
]

TARGET_STATES = [
    TargetState(
        property_types=[
            ("Density", 1),
        ],
        states=[
            State(
                temperature=298.15,
                pressure=101.325,
                mole_fractions=(1.0,),
            ),
        ],
    ),
    TargetState(
        property_types=[
            ("Density", 2),
            ("EnthalpyOfMixing", 2),
        ],
        states=[
            State(
                temperature=298.15,
                pressure=101.325,
                mole_fractions=(0.25, 0.75),
            ),
            State(
                temperature=298.15,
                pressure=101.325,
                mole_fractions=(0.5, 0.5),
            ),
            State(
                temperature=298.15,
                pressure=101.325,
                mole_fractions=(0.75, 0.25),
            ),
        ],
    ),
]

N_PROCESSES = 4


def curate_data_set(
    input_data_frame,
    property_type_filter: filtering.FilterByPropertyTypesSchema,
    allowed_elements: List[str],
) -> pandas.DataFrame:

    curation_schema = CurationWorkflowSchema(
        component_schemas=[
            # Filter out any measurements made for systems with more than
            # two components
            filtering.FilterByNComponentsSchema(n_components=[1, 2]),
            # Remove any duplicate data.
            filtering.FilterDuplicatesSchema(),
            # Filter out data points measured away from ambient conditions.
            filtering.FilterByTemperatureSchema(
                minimum_temperature=288.15, maximum_temperature=318.15
            ),
            filtering.FilterByPressureSchema(
                minimum_pressure=99.9, maximum_pressure=101.4
            ),
            # Retain only density and enthalpy of mixing data points which
            # have been measured for the same systems.
            property_type_filter,
            # Filter out long chain molecules (slower to simulate / converge) and 1, 3
            # carbonyl compounds where one of the carbonyls is a ketone (cases where
            # the enol form may be present in non-negligible amounts).
            filtering.FilterBySmirksSchema(
                smirks_to_exclude=[
                    # Long chain alkane /ether
                    "-".join(["[#6X4,#8X2]"] * 10),
                    # 1, 3 carbonyls with at least one ketone carbonyl.
                    "[#6](=[#8])-[#6](-[#1])(-[#1])-[#6](=[#8])-[#6]",
                ],
            ),
            # Filter out problematic molecules
            filtering.FilterBySmilesSchema(
                smiles_to_exclude=[
                    # Heavy water.
                    "[2H]O[2H]",
                    # Molecules which OpenMM misinterprets
                    "N[C@@H](CS)C(=O)O",
                    "CSCC[C@H](N)C(=O)O",
                    # Molecules which cause NaNs during simulations
                    "O=S(=O)(O)CCCN1CCOCC1",
                ]
            ),
            # Filter out systems where one component is in a significant excess.
            filtering.FilterByMoleFractionSchema(
                mole_fraction_ranges={2: [[(0.05, 0.95)]]}
            ),
            # Filter out any racemic mixtures
            filtering.FilterByRacemicSchema(),
            # Remove any substances measured for systems with undefined
            # stereochemistry
            filtering.FilterByStereochemistrySchema(),
            # Remove any measurements made for systems where any of the components
            # are charged.
            filtering.FilterByChargedSchema(),
            # Remove measurements made for ionic liquids
            filtering.FilterByIonicLiquidSchema(),
            # Remove any molecules containing elements that aren't currently of interest
            filtering.FilterByElementsSchema(allowed_elements=allowed_elements),
            # Retain only measurements made for substances which contain environments
            # of interest.
            filtering.FilterByEnvironmentsSchema(environments=CHEMICAL_ENVIRONMENTS),
            # Attempt to select a reasonable number of diverse substances
            selection.SelectSubstancesSchema(
                target_environments=CHEMICAL_ENVIRONMENTS,
                n_per_environment=5,
                per_property=False,
            ),
            # Select the data points for different compositions.
            selection.SelectDataPointsSchema(target_states=TARGET_STATES),
            # Filter out the density of water.
            filtering.FilterBySubstancesSchema(substances_to_exclude=[("O",)]),
        ]
    )

    return CurationWorkflow.apply(input_data_frame, curation_schema, N_PROCESSES)


def main():

    thermoml_data_frame = thermoml.ImportThermoMLData.apply(
        pandas.DataFrame(),
        thermoml.ImportThermoMLDataSchema(cache_file_name="thermoml.csv"),
        N_PROCESSES,
    )

    training_set_frame = curate_data_set(
        thermoml_data_frame,
        filtering.FilterByPropertyTypesSchema(
            property_types=[
                "Density",
                "EnthalpyOfMixing",
            ],
            n_components={
                "Density": [1, 2],
                "EnthalpyOfMixing": [2],
            },
            strict=True,
        ),
        ["C", "O", "N", "Cl", "Br", "H"],
    )

    training_set = DataSet.from_pandas(
        data_frame=training_set_frame,
        identifier="sage-train-v1",
        description="The main training set used for the Sage vdW refits.",
        authors=AUTHORS,
    )

    os.makedirs("data-sets", exist_ok=True)

    training_set.to_pandas().to_csv(
        os.path.join("data-sets", f"{training_set.id}.csv"), index=False
    )
    training_set.to_file(os.path.join("data-sets", f"{training_set.id}.json"))


if __name__ == "__main__":
    main()
