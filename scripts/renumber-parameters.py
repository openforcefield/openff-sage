"""A script that will update the id of each parameters in a force field such that
the id is of the form {prefix}{index} where {prefix} is the unique character associated
with a parameter handler and {index} is a sequential index starting from 1.

Certain parameters (i.e. those coming from the TIP3P water model) will be assigned
more descriptive id's so they can be more easily identified.
"""
import click
from openff.toolkit.typing.engines.smirnoff import ForceField


@click.option("-i", "--input", "input_file_name")
@click.option("-o", "--output", "output_file_name")
@click.command()
def main(input_file_name: str, output_file_name: str):

    force_field = ForceField(input_file_name)

    parameter_id_prefixes = {
        "Constraints": "c",
        "Bonds": "b",
        "Angles": "a",
        "ProperTorsions": "t",
        "ImproperTorsions": "i",
        "vdW": "n",
        "Electrostatics": "",
        "LibraryCharges": "",
        "ToolkitAM1BCC": ""
    }
    special_parameter_ids = {
        "Constraints": {
            "[#1:1]-[#8X2H2+0:2]-[#1]": "c-tip3p-H-O",
            "[#1:1]-[#8X2H2+0]-[#1:2]": "c-tip3p-H-O-H",
        },
        "LibraryCharges": {
            "[#1]-[#8X2H2+0:1]-[#1]": "q-tip3p-O",
            "[#1:1]-[#8X2H2+0]-[#1]": "q-tip3p-H",
        },
        "vdW": {
            "[#1]-[#8X2H2+0:1]-[#1]": "n-tip3p-O",
            "[#1:1]-[#8X2H2+0]-[#1]": "n-tip3p-H",
        },
    }

    for parameter_handler_name in force_field.registered_parameter_handlers:

        if parameter_handler_name not in parameter_id_prefixes:

            raise RuntimeError(
                f"The force field contains an unsupported handler: "
                f"{parameter_handler_name}."
            )

        parameter_handler = force_field[parameter_handler_name]

        for index, parameter in enumerate(parameter_handler.parameters):

            if (
                parameter_handler_name in special_parameter_ids
                and parameter.smirks in special_parameter_ids[parameter_handler_name]
            ):

                parameter.id = special_parameter_ids[parameter_handler_name][parameter.smirks]
                continue

            parameter_id_prefix = parameter_id_prefixes[parameter_handler_name]

            if parameter_id_prefix == "":
                continue

            parameter.id = f"{parameter_id_prefix}{index + 1}"

    force_field.to_file(output_file_name, "XML")

    # Sanity check that the force field can still be loaded.
    ForceField(output_file_name)


if __name__ == '__main__':
    main()
