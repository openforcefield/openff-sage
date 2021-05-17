# Physical Property Benchmark Data Sets

Within this directories are the scripts used to curate the experimental solvation / transfer free energy physical 
property data sets that the vdW force field parameters were tested against.

The data is sourced from both the FreeSolv and the Minnesota Solvation (MNSol) Databases. Due to how the MNSol is
licensed we do not include it here, however it can easily be obtained from the [main website](https://conservancy.umn.edu/handle/11299/213300).

## Curating the data sets

After obtaining a copy of the MNSol set:

1) Run the `mnsol-to-csv.py` script to convert the MNSol data from an XLS file to a more easily 
   manipulable CSV file.
   
2) Select a diverse and representative subset of the full MNSol by running the `curate-mnsol-test-set.py`.
   This will create both a `sage-mnsol-test-v1.json` and a `sage-mnsol-test-v1.csv` file in the `data-sets`
   directory which contain the selected test set data in two convenient formats. For completeness, we have 
   included a `sage-mnsol-test-v1.txt` file which contains exactly the substances used in our test set but 
   which excludes the actual values due to licensing concerns.
   
3) Select a subset of FreeSolv whereby each selected solute also appears as a solute in the selected MNSol
   set. This allows non-aques to aqueous transfer free energies to be computed. This will create both a 
   `sage-fsolv-test-v1.json` and a `sage-fsolv-test-v1.csv` file in the `data-sets` directory.
