#!/bin/bash

##SBATCH --job-name="full_parameter_sweep_test"
##SBATCH --time=1:00:00

# Initialize and Load Modules
# no modules to load when using MATLAB

echo "My task ID: " $LLSUB_RANK
echo "Number of Tasks: " $LLSUB_SIZE

mkdir data

matlab -batch "MyTaskID=$LLSUB_RANK;NumberOfTasks=$LLSUB_SIZE;nonlin_lookup"