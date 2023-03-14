Directories
 - fsolv*: Benchmark on FreeSolv database
 - mnsol*: Benchmark on Minnesota Solvation database
 
Each directory contains
 - benchmark.json: input file for the program nonbonded
 - estimation-options.json: input file for the program nonbonded, which would specify the protocol for calculating solvation free energies
 - forcefield.offxml: forcefield for which the properties are estimated
 - server-config.json: input file for the program nonbonded, which would specify the parallel worker configuration
 - test-set-collection.json: the reference data collected from the datasets mentioned, along with the source DOIs, std_errors, and other information
