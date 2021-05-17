# Physical Property Optimization Data Sets

Within this directories are the scripts used to curate the experimental density and enthalpy of mixing physical 
property data sets that the vdW force field parameters were trained against.

The data is sourced from the [ThermoML archive](https://www.nist.gov/mml/acmd/trc/thermoml) using the 
`openff-evaluator` packages data curation tools.

## Curating the data sets

1) Run the `curate-training-set.py` script to curate the training set. This will create both a 
   `sage-train-v1.csv` and a `sage-train-v1.json` file in the `data-sets` directory. The `sage-train-v1.json`
   file contains a JSON serialized `nonbonded.library.models.datasets.DataSet` object which can readily
   be deserialized by running `DataSet.parse_file("sage-train-v1.json")` from a python script.
