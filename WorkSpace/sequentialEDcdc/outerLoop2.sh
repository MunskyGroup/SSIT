#!/bin/bash
#SBATCH -n 16
#SBATCH --partition=cpu

module load matlab-2022b  # Load MATLAB module if necessary

# Loop over combinations of i, j, and k
for i in {3..3}; do
    for j in {1..2}; do
        for k in {1..1}; do
			LOG_FILE="log_${i}_${j}_${k}.txt"
			ERR_FILE="err_${i}_${j}_${k}.err"
        	srun matlab -nojvm -nodisplay -r "sequentialExptDesignBatchRunner($i, $j, $k); exit;" > "$LOG_FILE" 2> "$ERR_FILE" &
        done
    done
done

wait