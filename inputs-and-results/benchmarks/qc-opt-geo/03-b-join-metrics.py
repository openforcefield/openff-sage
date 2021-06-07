from typing import List

import click
import pandas


@click.command()
@click.argument(
    "file_names",
    nargs=-1,
    type=click.Path(exists=True, dir_okay=False, file_okay=True),
)
@click.option(
    "-o",
    "--output-path",
    "output_path",
    type=click.Path(exists=False, dir_okay=False),
)
def main(file_names: List[str], output_path: str):

    data_frame = pandas.concat(
        [pandas.read_csv(file_name) for file_name in file_names],
        ignore_index=True,
        sort=False,
    )

    data_frame = data_frame.sort_values(by=["SMILES", "Conformer Idx", "Force Field"])
    data_frame.to_csv(output_path, index=False)


if __name__ == "__main__":
    main()
