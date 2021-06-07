import os

import click
from openeye import oechem
from tqdm import tqdm


@click.command()
@click.option(
    "-i",
    "--input",
    "input_path",
    type=click.Path(exists=True, dir_okay=False),
    default="01-processed-qm.sdf",
)
@click.option(
    "-o", "--output-name", "output_name", type=click.STRING, default="01-processed-qm"
)
@click.option(
    "-dir",
    "--output-dir",
    "output_directory",
    type=click.Path(exists=False, dir_okay=True, file_okay=False),
    default="02-chunks",
)
@click.option("-n", "--n-chunks", "n_chunks", type=click.INT, default=250)
def main(input_path: str, output_name: str, output_directory: str, n_chunks: int):

    os.makedirs(output_directory, exist_ok=True)

    input_stream = oechem.oemolistream(input_path)

    n_molecules = sum(1 for _ in input_stream.GetOEGraphMols())
    input_stream.rewind()

    chunk_size = (n_molecules + n_chunks - 1) // n_chunks
    chunk_index = 1

    print(f"Chunking file N mols={n_molecules} N chunk={n_chunks} size={chunk_size}")

    output_stream = oechem.oemolostream(
        os.path.join(output_directory, f"{output_name}-{chunk_index}.sdf")
    )

    chunk_counter = 0

    for oe_molecule in tqdm(input_stream.GetOEGraphMols(), total=n_molecules):

        oe_molecule = oechem.OEMol(oe_molecule)

        if chunk_counter >= chunk_size:

            chunk_index += 1

            output_stream.close()

            output_stream.open(
                os.path.join(output_directory, f"{output_name}-{chunk_index}.sdf")
            )

            chunk_counter = 0

        chunk_counter += 1

        oechem.OEWriteMolecule(output_stream, oe_molecule)

    input_stream.close()
    output_stream.close()


if __name__ == "__main__":
    main()
