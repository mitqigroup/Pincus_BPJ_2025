#!/bin/bash

# Initialize and Load Modules
# no modules to load when using MATLAB

echo "My task ID: " $LLSUB_RANK
echo "Number of Tasks: " $LLSUB_SIZE

mkdir data

matlab -batch "MyTaskID=$LLSUB_RANK;NumberOfTasks=$LLSUB_SIZE;Toroid_Linear_Nonlin_param_sweep_supercloud"

