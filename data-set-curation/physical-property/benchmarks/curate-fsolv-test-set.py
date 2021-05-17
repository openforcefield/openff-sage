import os

import pandas
from nonbonded.library.models.authors import Author
from nonbonded.library.models.datasets import DataSet
from openff.evaluator.datasets.curation.components import filtering
from openff.evaluator.datasets.curation.components.freesolv import ImportFreeSolvSchema
from openff.evaluator.datasets.curation.workflow import (
    CurationWorkflow,
    CurationWorkflowSchema,
)
from openff.evaluator.utils.checkmol import ChemicalEnvironment

AUTHORS = [
    Author(
        name="Simon Boothroyd",
        email="simon.boothroyd@hotmail.com",
        institute="Boothroyd Scientific Consulting",
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

N_PROCESSES = 4


def main():

    # Extract the solutes present in the MNSol set.
    mnsol_data = pandas.read_csv(os.path.join("data-sets", "sage-mnsol-test-v1.csv"))

    solutes = [
        row["Component 1"] if row["Role 1"] == "Solute" else row["Component 2"]
        for _, row in mnsol_data.iterrows()
    ]

    print(f"Expected solutes ({len(solutes)}) - {solutes}")

    curation_schema = CurationWorkflowSchema(
        component_schemas=[
            ImportFreeSolvSchema(),
            # Retain on solutes prsent in the MNSol set
            filtering.FilterBySubstancesSchema(
                substances_to_include=[
                    (solute_smiles, "O") for solute_smiles in solutes
                ]
            ),
            filtering.FilterByStereochemistrySchema(),
        ]
    )
    # Apply the curation schema to yield the test set.
    full_test_set_frame = CurationWorkflow.apply(pandas.DataFrame(), curation_schema)

    print(f"Found solutes ({len(full_test_set_frame)})")

    test_set = DataSet.from_pandas(
        data_frame=full_test_set_frame,
        identifier="sage-fsolv-test-v1",
        description="A diverse set of aqueous solvation free energies to test the "
        "Sage line of force fields against. The solutes in this set where selected to "
        "match those in the `sage-mnsol-test-v1`.",
        authors=AUTHORS,
    )

    for i, entry in enumerate(test_set.entries):
        entry.id = i + 1

    os.makedirs("data-sets", exist_ok=True)

    test_set.to_pandas().to_csv(
        os.path.join("data-sets", f"{test_set.id}.csv"), index=False
    )
    test_set.to_file(os.path.join("data-sets", f"{test_set.id}.json"))


if __name__ == "__main__":
    main()
