#!/bin/bash

#SBATCH --job-name=calvariogram
#SBATCH --output=calvariogram.out
#SBATCH --ntasks=52
#SBATCH --partition=broadwl

# Load the default OpenMPI module.
module load openmpi

# Run the hellompi program with mpirun. The -n flag is not required;
# mpirun will automatically figure out the best configuration from the
# Slurm environment variables.
mpirun ./calvariogram
