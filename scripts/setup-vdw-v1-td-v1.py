import json
import os.path

from nonbonded.library.models.results import OptimizationResult
from openff.bespokefit.optimizers.forcebalance import ForceBalanceInputFactory
from openff.bespokefit.schema.fitting import OptimizationSchema
from openff.bespokefit.schema.optimizers import ForceBalanceSchema
from openff.bespokefit.schema.smirks import ProperTorsionSmirks
from openff.bespokefit.schema.smirnoff import ProperTorsionSettings
from openff.bespokefit.schema.targets import TorsionProfileTargetSchema
from openff.qcsubmit.results import TorsionDriveResultCollection


def main():

    torsion_training_set = TorsionDriveResultCollection.parse_file(
        "../data-set-curation/quantum-chemical/data-sets/1-2-0-td-set.json"
    )

    # Retrieve the FF with the fit vdW parameters and remove constraints as FB QM
    # targets do not support SMIRNOFF force fields which contain these.
    optimization_results = OptimizationResult.from_rest(
        project_id="openff-force-fields", study_id="sage", model_id="vdw-v1"
    )

    initial_force_field = optimization_results.refit_force_field.to_openff()
    initial_force_field.deregister_parameter_handler("Constraints")
    initial_force_field.to_file("vdw-v1.offxml")

    # Define the parameters to train
    with open(
        "../data-set-curation/quantum-chemical/data-sets/"
        "1-2-0-td-set-torsion-smirks.json"
    ) as file:

        torsion_smirks = json.load(file)

    torsion_handler = initial_force_field.get_parameter_handler("ProperTorsions")

    target_parameters = [
        ProperTorsionSmirks(
            smirks=smirks,
            attributes={
                f"k{i + 1}" for i in range(len(torsion_handler.parameters[smirks].k))
            },
        )
        for smirks in torsion_smirks["ProperTorsions"]
    ]

    # Define the full schema for the optimization.
    optimization_schema = OptimizationSchema(
        id="vdw-v1-td-v1",
        initial_force_field=os.path.abspath("vdw-v1.offxml"),
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
        ],
        # Define the parameters to refit and the priors to place on them.
        target_parameters=target_parameters,
        parameter_settings=[ProperTorsionSettings(target="k", prior=2.0)],
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

    # Remove the temporary force field
    os.unlink("vdw-v1.offxml")


if __name__ == "__main__":
    main()
