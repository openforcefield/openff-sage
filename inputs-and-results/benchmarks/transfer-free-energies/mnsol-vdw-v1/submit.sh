#!/bin/bash
#
# Set the job name and wall time / memory limit
#BSUB -J bench
#BSUB -W 168:00
#BSUB -M 8
#
# Set the output and error output paths.
#BSUB -o  %J.o
#BSUB -e  %J.e
#
#BSUB -q cpuqueue


. ~/.bashrc

# Use the right conda environment
conda activate openff-force-fields
conda env export > conda_env.yaml

# Run the commands
nonbonded benchmark run --restart true
nonbonded benchmark analyze
