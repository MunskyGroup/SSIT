#!/bin/bash

# Loop over combinations of i, j, and k
for i in {8..8}; do
    for j in {1..2}; do
        for k in {2..2}; do
            echo "running job: " "$i" "$j" "$k"
			sbatch batch_run_sequential_experiment_analyses.sh "$i" "$j" "$k"
        done
    done
done