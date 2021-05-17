import copy
import functools
import json
import logging
from collections import defaultdict
from multiprocessing import Pool
from tempfile import NamedTemporaryFile

from openff.qcsubmit.results import TorsionDriveResultCollection
from openff.qcsubmit.results.filters import (
    ConnectivityFilter,
    HydrogenBondFilter,
    RecordStatusFilter,
    ResultRecordFilter,
    SMARTSFilter,
)
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils import UndefinedStereochemistryError
from qcportal import FractalClient
from qcportal.models.records import RecordStatusEnum
from tqdm import tqdm


class UndefinedStereoFilter(ResultRecordFilter):
    def _filter_function(self, result, record, molecule) -> bool:

        has_stereochemistry = True

        molecule = copy.deepcopy(molecule)
        molecule._conformers = [molecule.conformers[0]]

        try:

            with NamedTemporaryFile(suffix=".sdf") as file:
                molecule.to_file(file.name, "SDF")
                molecule.from_file(file.name)

        except UndefinedStereochemistryError:
            has_stereochemistry = False

        return has_stereochemistry


def label_ids(record_and_molecule, force_field):

    record, molecule = record_and_molecule

    full_labels = force_field.label_molecules(molecule.to_topology())[0]
    torsion_labels = full_labels["ProperTorsions"]

    parameter_ids = set()

    for indices, parameter in torsion_labels.items():

        if {*indices[1:3]} != {*record.keywords.dihedrals[0][1:3]}:
            continue

        parameter_ids.add(parameter.id)

    return [*parameter_ids]


def main():

    logging.basicConfig(level=logging.INFO)

    # Pull down the main amide torsion drive set and filter out any records which
    # have not completed or which inadvertently contain intra-molecular h-bonds.
    results = TorsionDriveResultCollection.from_server(
        client=FractalClient(),
        datasets=[
            "OpenFF Gen3 Torsion Set v1.0",
            "OpenFF Amide Torsion Set v1.0",
        ],
        spec_name="default",
    )

    results = results.filter(
        SMARTSFilter(
            # Filter out unusual chemistries.
            smarts_to_exclude=[
                "[#7]=[#7X3]-[#8X1H0]",
                "[#7]-[#7]=[#8]",
                "[#15]-[#7](~[#8])(~[#8])",
            ]
        ),
        RecordStatusFilter(status=RecordStatusEnum.complete),
        HydrogenBondFilter(method="baker-hubbard"),
        ConnectivityFilter(tolerance=1.2),
        UndefinedStereoFilter(),
    )

    with open("data-sets/2-0-0-td-set-v2.json", "w") as file:
        file.write(results.json())

    results = TorsionDriveResultCollection.parse_file("data-sets/2-0-0-td-set-v2.json")

    # Print out coverage information.
    force_field = ForceField("openff-1.3.0.offxml")

    coverage = defaultdict(int)

    with Pool(16) as pool:

        for parameter_ids in tqdm(
            pool.imap(
                functools.partial(label_ids, force_field=force_field),
                results.to_records(),
            ),
            total=results.n_results,
        ):

            for parameter_id in parameter_ids:
                coverage[parameter_id] += 1

    print(coverage)

    # Save out the SMIRKS which should be trained against this set.
    with open("data-sets/2-0-0-td-set-v2-torsion-smirks.json", "w") as file:

        json.dump(
            [
                force_field.get_parameter_handler("ProperTorsions")
                .get_parameter({"id": parameter_id})[0]
                .smirks
                for parameter_id, count in coverage.items()
                if count >= 5
            ],
            file,
        )


if __name__ == "__main__":
    main()
