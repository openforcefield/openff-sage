import os
from glob import glob

import click
from openeye import oechem
from tqdm import tqdm


@click.command()
@click.option(
    "-i",
    "--input-name",
    "input_name",
    type=click.STRING,
)
@click.option(
    "-dir",
    "--input-dir",
    "input_directory",
    type=click.Path(exists=True, dir_okay=True, file_okay=False),
    default="02-outputs",
)
@click.option(
    "-o",
    "--output-path",
    "output_path",
    type=click.Path(exists=False, dir_okay=False),
)
def main(input_name: str, input_directory: str, output_path: str):

    output_stream = oechem.oemolostream(output_path)

    n_inputs = len(glob(os.path.join(input_directory, f"{input_name}-*.sdf")))

    for i in tqdm(range(1, n_inputs + 1)):

        input_stream = oechem.oemolistream(
            os.path.join(input_directory, f"{input_name}-{i}.sdf")
        )

        for oe_molecule in input_stream.GetOEGraphMols():

            oe_molecule = oechem.OEGraphMol(oe_molecule)
            oechem.OEWriteMolecule(output_stream, oe_molecule)

        input_stream.close()

    output_stream.close()


if __name__ == "__main__":
    main()
