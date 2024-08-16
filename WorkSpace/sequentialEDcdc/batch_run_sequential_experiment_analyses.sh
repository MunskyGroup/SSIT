#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --partition=cpu

# Assign input arguments to variables
i=$1
j=$2
k=$3

module load matlab-2022b  # Load MATLAB module if necessary

# Run MATLAB job with input arguments
LOG_FILE="log_${i}_${j}_${k}.log"
ERR_FILE="err_${i}_${j}_${k}.err"
srun matlab -nojvm -nodisplay -r "sequentialExptDesignBatchRunner($i, $j, $k); exit;" > "$LOG_FILE" 2> "$ERR_FILE" &
wait