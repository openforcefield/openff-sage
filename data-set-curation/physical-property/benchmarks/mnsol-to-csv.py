import json

import pandas
from openff.evaluator import unit
from openff.evaluator.datasets import (
    MeasurementSource,
    PhysicalPropertyDataSet,
    PropertyPhase,
)
from openff.evaluator.properties import SolvationFreeEnergy
from openff.evaluator.substances import Component, ExactAmount, MoleFraction, Substance
from openff.evaluator.thermodynamics import ThermodynamicState


def main():

    original_data = pandas.read_excel("MNSol_alldata.xls")

    # Load in a map between most of the common names found in the MNSol set and the
    # corresponding SMILES patterns generated using pubchem and NIH cactus.
    with open("name-to-smiles.json") as file:
        name_to_smiles = json.load(file)

    solute_smiles = original_data["SoluteName"].apply(
        lambda name: name_to_smiles.get(name, None)
    )
    solvent_smiles = original_data["Solvent"].apply(
        lambda name: name_to_smiles.get(name, None)
    )

    original_data["Solute SMILES"] = solute_smiles
    original_data["Solvent SMILES"] = solvent_smiles

    original_data = original_data[
        original_data["Solvent SMILES"].notna() & original_data["Solute SMILES"].notna()
    ]

    properties = []

    for index, data_row in original_data.iterrows():

        substance = Substance()
        substance.add_component(
            Component(smiles=data_row["Solvent SMILES"], role=Component.Role.Solvent),
            MoleFraction(1.0),
        )
        substance.add_component(
            Component(smiles=data_row["Solute SMILES"], role=Component.Role.Solute),
            ExactAmount(1),
        )

        solvation_free_energy = SolvationFreeEnergy(
            thermodynamic_state=ThermodynamicState(
                temperature=298.0 * unit.kelvin, pressure=1.0 * unit.atmosphere
            ),
            phase=PropertyPhase.Liquid,
            substance=substance,
            value=data_row["DeltaGsolv"] * unit.kilocalorie / unit.mole,
            uncertainty=0.2 * unit.kilocalorie / unit.mole,
            source=MeasurementSource(doi="10.13020/3eks-j059"),
        )

        properties.append(solvation_free_energy)

    openff_data_set = PhysicalPropertyDataSet()
    openff_data_set.add_properties(*properties, validate=False)
    openff_data_set.to_pandas().to_csv("off_converted.csv", index=False)


if __name__ == "__main__":
    main()
