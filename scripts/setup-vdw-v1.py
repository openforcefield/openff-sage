import os

from nonbonded.library.models.engines import ForceBalance
from nonbonded.library.models.forcefield import ForceField, Parameter
from nonbonded.library.models.projects import Optimization
from nonbonded.library.models.targets import EvaluatorTarget
from nonbonded.library.utilities.environments import ChemicalEnvironment
from openff.toolkit.typing.engines import smirnoff

PROJECT_ID = "openff-force-fields"
STUDY_ID = "sage"

ANALYSIS_ENVIRONMENTS = [
    ChemicalEnvironment.Aqueous,
    ChemicalEnvironment.SecondaryAmine,
    ChemicalEnvironment.CarboxylicAcidSecondaryAmide,
    ChemicalEnvironment.AlkylBromide,
    ChemicalEnvironment.Alcohol,
    ChemicalEnvironment.Aromatic,
    ChemicalEnvironment.ArylChloride,
    ChemicalEnvironment.Thiol,
    ChemicalEnvironment.CarboxylicAcidTertiaryAmide,
    ChemicalEnvironment.AlkylChloride,
    ChemicalEnvironment.CarboxylicAcidEster,
    ChemicalEnvironment.Thiourea,
    ChemicalEnvironment.Acetal,
    ChemicalEnvironment.Thioether,
    ChemicalEnvironment.TertiaryAmine,
    ChemicalEnvironment.Ketone,
    ChemicalEnvironment.Disulfide,
    ChemicalEnvironment.Aldehyde,
    ChemicalEnvironment.PrimaryAmine,
    ChemicalEnvironment.Heterocycle,
    ChemicalEnvironment.Ether,
    ChemicalEnvironment.Alkene,
    ChemicalEnvironment.Alkane,
    ChemicalEnvironment.Sulfone,
]


def create_tip3p_force_field() -> str:

    return (
        '<SMIRNOFF version="0.3" aromaticity_model="OEAroModel_MDL">'
        "  <Date>2020-09-01</Date>"
        "  <Author>S. Boothroyd</Author>"
        '  <Constraints version="0.3">'
        '    <Constraint smirks="[#1:1]-[#8X2H2+0:2]-[#1]" id="c1" distance="0.9572 * angstrom"/>'
        '    <Constraint smirks="[#1:1]-[#8X2H2+0]-[#1:2]" id="c2" distance="1.5139006545247014 * angstrom"/>'
        "  </Constraints>"
        '  <vdW version="0.3" potential="Lennard-Jones-12-6" combining_rules="Lorentz-Berthelot" scale12="0.0" scale13="0.0" scale14="0.5" scale15="1.0" cutoff="9.0 * angstrom" switch_width="1.0 * angstrom" method="cutoff">'
        '    <Atom smirks="[#1]-[#8X2H2+0:1]-[#1]" epsilon="0.1521 * mole**-1 * kilocalorie" id="n1" sigma="3.1507 * angstrom"/>'
        '    <Atom smirks="[#1:1]-[#8X2H2+0]-[#1]" epsilon="0 * mole**-1 * kilocalorie" id="n2" sigma="1 * angstrom"/>'
        "  </vdW>"
        '  <LibraryCharges version="0.3">'
        '    <LibraryCharge smirks="[#1]-[#8X2H2+0:1]-[#1]" charge1="-0.834 * elementary_charge"/>'
        '    <LibraryCharge smirks="[#1:1]-[#8X2H2+0]-[#1]" charge1="0.417 * elementary_charge"/>'
        "  </LibraryCharges>"
        "</SMIRNOFF>"
    )


def main():

    optimization = Optimization(
        project_id=PROJECT_ID,
        study_id=STUDY_ID,
        id="vdw-v1",
        name="VdW Parameters V1",
        description=(
            "An optimization of the vdW parameters of the `openff-1.3.0` force "
            "field against a training set of physical property data."
        ),
        engine=ForceBalance(
            priors={
                "vdW/Atom/epsilon": 0.1,
                "vdW/Atom/rmin_half": 1.0,
            }
        ),
        targets=[
            EvaluatorTarget(
                id="phys-prop",
                denominators={
                    "Density": "0.05 g / ml",
                    "EnthalpyOfMixing": "1.6 kJ / mol",
                },
                data_set_ids=["sage-train-v1"],
            )
        ],
        force_field=ForceField.from_openff(
            smirnoff.ForceField("openff-1.3.0.offxml", create_tip3p_force_field())
        ),
        parameters_to_train=[
            Parameter(handler_type="vdW", attribute_name=attribute, smirks=smirks)
            for attribute in ["epsilon", "rmin_half"]
            for smirks in [
                "[#16:1]",
                "[#17:1]",
                "[#1:1]-[#6X3]",
                "[#1:1]-[#6X3](~[#7,#8,#9,#16,#17,#35])~[#7,#8,#9,#16,#17,#35]",
                "[#1:1]-[#6X3]~[#7,#8,#9,#16,#17,#35]",
                "[#1:1]-[#6X4]",
                "[#1:1]-[#6X4]-[#7,#8,#9,#16,#17,#35]",
                "[#1:1]-[#7]",
                "[#1:1]-[#8]",
                "[#35:1]",
                "[#6:1]",
                "[#6X4:1]",
                "[#7:1]",
                "[#8:1]",
                "[#8X2H0+0:1]",
                "[#8X2H1+0:1]",
            ]
        ],
        analysis_environments=ANALYSIS_ENVIRONMENTS,
        max_iterations=15,
    )

    with open(
        os.path.join(os.pardir, "schemas", "optimizations", f"{optimization.id}.json"),
        "w",
    ) as file:
        file.write(optimization.json())

    optimization.upload()


if __name__ == "__main__":
    main()
