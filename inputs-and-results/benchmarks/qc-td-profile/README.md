# QC Torsion Profile Benchmarks

This directory contains the scripts needed to compare the QM and MM torsion profiles for a set
of force fields against a QCArchive torsion drive collection (e.g. `OpenFF-benchmark-ligand-fragments-v1.0`).

## Running the benchmarks

1) Run the `01-setup.py` script to retrive the QC torsion drive data and compute the associated MM
   torsion profiles and residuals, e.g.:
   
   ```shell
   python 01-setup.py --data-set "OpenFF-benchmark-ligand-fragments-v1.0"
   ```
   
2) Plot the metrics about the torsion drives, and well as the torsion profiles themselves by running the 
   `02-plot.py` script.
