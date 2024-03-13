#!/bin/bash
#SBATCH --job-name=my_matlab_job
#SBATCH --output=my_matlab_job.out
#SBATCH --error=my_matlab_job.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=all

module load matlab-2022b  # Load MATLAB module if necessary

# Loop over combinations of i and j
for i in {1..3}; do
    for j in {1..2}; do
    	for k in {1..3}; do
    	    LOG_FILE="log_${i}_${j}_${k}.txt"
        	# Submit MATLAB job to Slurm for each combination
        	matlab -nojvm -nodisplay -r "sequentialExptDesignBatchRunner($i, $j, $k); exit;" > "$LOG_FILE" &
     done
done
