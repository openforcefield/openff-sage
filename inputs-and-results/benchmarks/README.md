# Benchmarks

Contained within this directory are the main inputs (or instructions on how they may be generated)
for performing benchmarks of the force fields produced during Sages' development against QC, protein-ligand 
binding free energy, and solvation / transfer free energy data.

## Manifest

* `qc-opt-geo` - a directory containing the inputs required to compute the ddE, RMSD, and TDF metrics of 
  each produced force field against data sets retrieved from the QCArchive.

* `qc-opt-geo` - a directory containing scripts to compute and compare a set of computed MM torsion profiles
  to a QC set retrieved from QCArchive.

* `protein-ligand` - a directory containing the inputs required to perform relative bind free energy
  calculations of the TYK2 JACS ligand set using the `perses` package.
  
* `transfer-free-energy` - a directory containing the inputs required to perform solvation / transfer free energy
  calculations of a diverse set of solutes / solvents using the `openff-evaluator` package.
