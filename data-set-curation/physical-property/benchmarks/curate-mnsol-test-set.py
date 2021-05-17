import os
from typing import Dict

import pandas
from nonbonded.library.models.authors import Author
from nonbonded.library.models.datasets import DataSet
from openff.evaluator.datasets.curation.components import filtering, selection
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


def available_solvents() -> Dict[str, str]:
    return {
        "2methylpyridine": "Cc1ccccn1",
        "4methyl2pentanone": "CC(C)CC(=O)C",
        "aceticacid": "CC(=O)O",
        "acetonitrile": "CC#N",
        "acetophenone": "CC(=O)c1ccccc1",
        "aniline": "c1ccc(cc1)N",
        "anisole": "COc1ccccc1",
        "benzene": "c1ccccc1",
        "benzonitrile": "c1ccc(cc1)C#N",
        "benzylalcohol": "c1ccc(cc1)CO",
        "bromobenzene": "c1ccc(cc1)Br",
        "bromoethane": "CCBr",
        "bromoform": "C(Br)(Br)Br",
        "bromooctane": "CCCCCCCCBr",
        "butanol": "CCCCO",
        "butanone": "CCC(=O)C",
        "butylacetate": "CCCCOC(=O)C",
        "butylbenzene": "CCCCc1ccccc1",
        "carbondisulfide": "C(=S)=S",
        "chlorobenzene": "c1ccc(cc1)Cl",
        "chloroform": "C(Cl)(Cl)Cl",
        "chlorohexane": "CCCCCCCl",
        "cyclohexane": "C1CCCCC1",
        "cyclohexanone": "C1CCC(=O)CC1",
        "decalin": "C1CCC2CCCCC2C1",
        "decane": "CCCCCCCCCC",
        "decanol": "CCCCCCCCCCO",
        "dibromoethane": "CC(Br)Br",
        "dibutylether": "CCCCOCCCC",
        "dichloroethane": "CC(Cl)Cl",
        "diethylether": "CCOCC",
        "diisopropylether": "CC(C)OC(C)C",
        "dimethylacetamide": "CC(=O)N(C)C",
        "dimethylformamide": "CN(C)C=O",
        "dimethylpyridine": "Cc1cccnc1C",
        "dimethylsulfoxide": "CS(=O)C",
        "dodecane": "CCCCCCCCCCCC",
        "ethanol": "CCO",
        "ethoxybenzene": "CCOc1ccccc1",
        "ethylacetate": "CCOC(=O)C",
        "ethylbenzene": "CCc1ccccc1",
        "fluorobenzene": "c1ccc(cc1)F",
        "heptane": "CCCCCCC",
        # "heptanol": "CCCCCCCO",
        "hexadecane": "CCCCCCCCCCCCCCCC",
        "hexadecyliodide": "CCCCCCCCCCCCCCCCI",
        "hexane": "CCCCCC",
        "hexanol": "CCCCCCO",
        "iodobenzene": "c1ccc(cc1)I",
        "isobutanol": "CC(C)CO",
        "isooctane": "CCCCCC(C)C",
        "isopropanol": "CC(C)O",
        "isopropylbenzene": "CC(C)c1ccccc1",
        "isopropyltoluene": "Cc1ccccc1C(C)C",
        "mesitylene": "Cc1cc(cc(c1)C)C",
        "methoxyethanol": "CC(O)OC",
        "methylenechloride": "C(Cl)Cl",
        "methylformamide": "CNC=O",
        "nitrobenzene": "c1ccc(cc1)[N+](=O)[O-]",
        "nitroethane": "CC[N+](=O)[O-]",
        "nitromethane": "C[N+](=O)[O-]",
        # "nonane": "CCCCCCCCC",
        # "nonanol": "CCCCCCCCCO",
        "octane": "CCCCCCCC",
        # "octanol": "CCCCCCCCO",
        # "pentadecane": "CCCCCCCCCCCCCCC",
        "pentane": "CCCCC",
        "pentanol": "CCCCCO",
        "perfluorobenzene": "c1(c(c(c(c(c1F)F)F)F)F)F",
        "phenylether": "c1ccc(cc1)Oc2ccccc2",
        "propanol": "CCCO",
        "pyridine": "c1ccncc1",
        "secbutanol": "CCC(C)O",
        "secbutylbenzene": "CCC(C)c1ccccc1",
        "tbutylbenzene": "CC(C)(C)C1=CC=CC=C1",
        "tetrachloroethene": "C(=C(Cl)Cl)(Cl)Cl",
        "tetrahydrofuran": "C1CCOC1",
        "tetrahydrothiophenedioxide": "C1CCS(=O)(=O)C1",
        "tetralin": "c1ccc2c(c1)CCCC2",
        "toluene": "Cc1ccccc1",
        "tributylphosphate": "CCCCOP(=O)(OCCCC)OCCCC",
        "triethylamine": "CCN(CC)CC",
        "trimethylbenzene": "Cc1cccc(c1C)C",
        "undecane": "CCCCCCCCCCC",
        # "water": "O",
        "xylene": "Cc1ccccc1C",
        # "methanol": "CO",
    }


def curate_experimental_set(solvent_smiles: str) -> pandas.DataFrame:
    """Curate the set of physical properties to use for the training set."""
    curation_schema = CurationWorkflowSchema(
        component_schemas=[
            # Remove measurements of the solvent with itself.
            filtering.FilterBySubstancesSchema(
                substances_to_exclude=[(solvent_smiles, solvent_smiles)]
            ),
            # Filter out long chain molecules (slower to simulate / converge) and 1, 3
            # carbonyl compounds where one of the carbonyls is a ketone (cases where
            # the enol form may be present in non-negligible amounts).
            filtering.FilterBySmirksSchema(
                smirks_to_exclude=[
                    # Long chain alkane /ether
                    "-".join(["[#6X4,#8X2]"] * 10),
                    # 1, 3 carbonyls with at least one ketone carbonyl.
                    "[#6](=[#8])-[#6](-[#1])(-[#1])-[#6](=[#8])-[#6]",
                    # Things which disociate
                    "[H]-[Cl]",
                    "[H]-[Br]",
                ],
            ),
            # Remove any substances measured for systems with undefined
            # stereochemistry
            filtering.FilterByStereochemistrySchema(),
            # Remove any measurements made for systems where any of the components
            # are charged.
            filtering.FilterByChargedSchema(),
            # Remove measurements made for ionic liquids
            filtering.FilterByIonicLiquidSchema(),
            # Remove any molecules containing unwanted elements
            filtering.FilterByElementsSchema(
                allowed_elements=["C", "O", "N", "Cl", "Br", "H", "S", "P"]
            ),
            # Retain only measurements made for substances which contain environments
            # of interest.
            filtering.FilterByEnvironmentsSchema(environments=CHEMICAL_ENVIRONMENTS),
            # Attempt to select a reasonable number of diverse substances
            selection.SelectSubstancesSchema(
                target_environments=CHEMICAL_ENVIRONMENTS,
                n_per_environment=2,
                per_property=False,
            ),
        ]
    )
    # Apply the curation schema to yield the test set.
    full_mnsol_data = pandas.read_csv("mnsol/off_converted.csv")

    solvent_data = full_mnsol_data[
        (
            (full_mnsol_data["Component 1"] == solvent_smiles)
            & full_mnsol_data["Mole Fraction 1"].notna()
        )
        | (
            (full_mnsol_data["Component 2"] == solvent_smiles)
            & full_mnsol_data["Mole Fraction 2"].notna()
        )
    ]

    try:
        return CurationWorkflow.apply(solvent_data, curation_schema, N_PROCESSES)
    except TypeError:

        return pandas.DataFrame()


def main():

    per_solvent_test_sets = {
        solvent_smiles: curate_experimental_set(solvent_smiles)
        for _, solvent_smiles in available_solvents().items()
    }

    full_test_set_frame = pandas.concat(
        per_solvent_test_sets.values(), sort=False, ignore_index=True
    )

    test_set = DataSet.from_pandas(
        data_frame=full_test_set_frame,
        identifier="sage-mnsol-test-v1",
        description="A diverse set of non-aqueous solvation free energies to test the "
        "Sage line of force fields against.",
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
