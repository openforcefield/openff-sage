import json
import os.path

from nonbonded.library.models.results import OptimizationResult
from openff.bespokefit.optimizers.forcebalance import ForceBalanceInputFactory
from openff.bespokefit.schema.fitting import OptimizationSchema
from openff.bespokefit.schema.optimizers import ForceBalanceSchema
from openff.bespokefit.schema.smirks import AngleSmirks, BondSmirks, ProperTorsionSmirks
from openff.bespokefit.schema.smirnoff import (
    AngleAngleSettings,
    AngleForceSettings,
    BondForceSettings,
    BondLengthSettings,
    ProperTorsionSettings,
)
from openff.bespokefit.schema.targets import (
    OptGeoTargetSchema,
    TorsionProfileTargetSchema,
)
from openff.qcsubmit.results import (
    OptimizationResultCollection,
    TorsionDriveResultCollection,
)
from openff.toolkit.typing.engines.smirnoff import ForceField


def main():

    torsion_training_set = TorsionDriveResultCollection.parse_file(
        "../data-set-curation/quantum-chemical/data-sets/1-2-0-td-set.json"
    )
    optimization_training_set = OptimizationResultCollection.parse_file(
        "../data-set-curation/quantum-chemical/data-sets/1-2-0-opt-set.json"
    )

    # Define the parameters to train
    original_force_field = ForceField("openff-1.3.0.offxml")
    initial_force_field = ForceField("openff-1-3-0-msm.offxml")

    for parameter in original_force_field["Bonds"].parameters:
        initial_force_field["Bonds"].get_parameter({"smirks": parameter.smirks})[
            0
        ].length = parameter.length
    for parameter in original_force_field["Angles"].parameters:
        initial_force_field["Angles"].get_parameter({"smirks": parameter.smirks})[
            0
        ].angle = parameter.angle

    initial_force_field.to_file("openff-1-3-0-msm.offxml")

    with open(
        "../data-set-curation/quantum-chemical/data-sets/"
        "1-2-0-opt-set-valence-smirks.json"
    ) as file:

        valence_smirks = json.load(file)

    with open(
        "../data-set-curation/quantum-chemical/data-sets/"
        "1-2-0-td-set-torsion-smirks.json"
    ) as file:

        torsion_smirks = json.load(file)

    target_parameters = [
        *[
            BondSmirks(smirks=smirks, attributes={"k", "length"})
            for smirks in valence_smirks["Bonds"]
        ],
        *[
            AngleSmirks(smirks=smirks, attributes={"k", "angle"})
            for smirks in valence_smirks["Angles"]
        ],
        *[
            ProperTorsionSmirks(
                smirks=smirks,
                attributes={
                    f"k{i + 1}"
                    for i in range(
                        len(
                            initial_force_field.get_parameter_handler("ProperTorsions")
                            .parameters[smirks]
                            .k
                        )
                    )
                },
            )
            for smirks in torsion_smirks["ProperTorsions"]
        ],
    ]

    # Define the full schema for the optimization.
    optimization_schema = OptimizationSchema(
        id="vdw-v1-ms-v1-td-opt-v1",
        initial_force_field=os.path.abspath("openff-1-3-0-msm.offxml"),
        # Define the optimizer / ForceBalance specific settings.
        optimizer=ForceBalanceSchema(
            max_iterations=50,
            step_convergence_threshold=0.01,
            objective_convergence_threshold=0.1,
            gradient_convergence_threshold=0.1,
            n_criteria=2,
            initial_trust_radius=-1.0,
            extras={"wq_port": "55125", "asynchronous": "True"},
        ),
        # Define the torsion profile targets to fit against.
        targets=[
            TorsionProfileTargetSchema(
                reference_data=torsion_training_set,
                energy_denominator=1.0,
                energy_cutoff=5.0,
                extras={"remote": "1"},
            ),
            OptGeoTargetSchema(
                reference_data=optimization_training_set,
                weight=0.1,
                extras={"batch_size": 30, "remote": "1"},
            ),
        ],
        # Define the parameters to refit and the priors to place on them.
        target_parameters=target_parameters,
        parameter_settings=[
            BondForceSettings(prior=1.0e01),
            BondLengthSettings(prior=1.0e-01),
            AngleForceSettings(prior=1.0e01),
            AngleAngleSettings(prior=2.0e01),
            ProperTorsionSettings(target="k", prior=1.0),
        ],
    )

    with open(
        os.path.join(
            os.pardir, "schemas", "optimizations", f"{optimization_schema.id}.json"
        ),
        "w",
    ) as file:
        file.write(optimization_schema.json())

    # Generate the ForceBalance inputs
    ForceBalanceInputFactory.generate(
        os.path.join(
            os.pardir, "inputs-and-results", "optimizations", optimization_schema.id
        ),
        optimization_schema,
    )


if __name__ == "__main__":
    main()
