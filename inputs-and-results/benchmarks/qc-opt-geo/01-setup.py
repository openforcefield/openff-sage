import logging
from collections import defaultdict
from multiprocessing import Pool
from typing import List, Optional

import click
import qcelemental
from openeye import oechem
from openff.qcsubmit.results import OptimizationResultCollection
from openff.qcsubmit.results.filters import (
    CMILESResultFilter,
    ConnectivityFilter,
    RecordStatusFilter,
    ResultRecordFilter,
)
from openff.toolkit.topology import Molecule
from openff.toolkit.typing.engines.smirnoff import ForceField
from openff.toolkit.utils import UndefinedStereochemistryError
from qcportal import FractalClient
from qcportal.models.records import RecordStatusEnum

N_PROCESSES = 24


class InvalidCMILESFilter(CMILESResultFilter):
    def _filter_function(self, result) -> bool:

        try:
            Molecule.from_mapped_smiles(result.cmiles, allow_undefined_stereo=True)
        except BaseException:
            logging.exception(f"skipping {result.record_id}")
            return False

        return True


def _process_molecule(record_and_molecule) -> Optional[oechem.OEMol]:
    """Convert a QC record and its associated molecule into an OE molecule which
    has been tagged with the associated SMILES, final energy and record id."""

    force_field = ForceField("openff-1.3.0.offxml")

    record, molecule = record_and_molecule

    print(f"processing {record.id}")

    oe_molecule = molecule.to_openeye()
    # Re-perceive the stereochemistry from the final conformer.
    oechem.OE3DToInternalStereo(oe_molecule)

    try:
        molecule = Molecule.from_openeye(oe_molecule)
    except UndefinedStereochemistryError:
        print(f"skipping {record.id} - un-perceivable stereo")
        return None

    try:
        force_field.create_openmm_system(molecule.to_topology())
    except BaseException:
        print(f"skipping {record.id} - cannot parameterize")
        return None

    final_energy = record.get_final_energy() * qcelemental.constants.hartree2kcalmol

    # add name and energy tag to the mol
    oechem.OESetSDData(oe_molecule, "SMILES QCArchive", molecule.to_smiles())
    oechem.OESetSDData(oe_molecule, "Energy QCArchive", str(final_energy))
    oechem.OESetSDData(oe_molecule, "Record QCArchive", str(record.id))

    return oe_molecule


@click.command()
@click.argument(
    "data_set",
    nargs=1,
    type=click.STRING,
)
def main(data_set):
    """Process and store optimized QC geometries from a QCArchive dataset.

    The DATA_SET should likely be one of:

    \b
    * "OpenFF Full Optimization Benchmark 1"
    * "OpenFF Industry Benchmark Season 1 v1.0"

    The processed molecules tagged with information from the QC record, including the
    CMILES and QC energy, will be stored in a new `01-processed-qm.sdf` file and
    additionally information about the included records will be stored in
    `01-processed-qm.json`.
    """

    print("1a) Parsing collection")

    client = FractalClient()

    result_collection = OptimizationResultCollection.from_server(
        client=client,
        datasets=data_set,
        spec_name="default",
    )

    print("1b) Filtering unfinished / undesirable results")

    result_collection = result_collection.filter(
        InvalidCMILESFilter(),
        RecordStatusFilter(status=RecordStatusEnum.complete),
        ConnectivityFilter(),
    )

    print("1c) Retrieving QC records")

    records_and_molecules = result_collection.to_records()

    print(f"1d) Grouping molecules (N={len(records_and_molecules)})")

    grouped_molecules = defaultdict(list)

    for record, molecule in records_and_molecules:

        molecule = molecule.canonical_order_atoms()

        smiles = molecule.to_smiles(isomeric=False, explicit_hydrogens=True)
        grouped_molecules[smiles].append((record, molecule))

    print(f"1e) Filtering unprocessable molecules (N={len(records_and_molecules)})")

    with Pool(processes=N_PROCESSES) as pool:

        processed_oe_molecules: List[oechem.OEMol] = [
            oe_molecule
            for oe_molecule in pool.map(
                _process_molecule,
                (pair for values in grouped_molecules.values() for pair in values),
            )
            if oe_molecule is not None
        ]

    print(f"1f) Writing processed molecules to sdf (N={len(processed_oe_molecules)})")

    output_steam = oechem.oemolostream("01-processed-qm.sdf")

    final_record_ids = set()

    for i, oe_molecule in enumerate(processed_oe_molecules):

        final_record_ids.add(oechem.OEGetSDData(oe_molecule, "Record QCArchive"))

        oe_molecule.SetTitle(f"full_{i + 1}")
        oechem.OEWriteConstMolecule(output_steam, oe_molecule)

    output_steam.close()

    result_collection.entries[client.address] = [
        entry
        for entry in result_collection.entries[client.address]
        if entry.record_id in final_record_ids
    ]

    with open("01-processed-qm.json", "w") as file:
        file.write(result_collection.json())


if __name__ == "__main__":
    main()
