#!/bin/bash
#!/bin/bash
#SBATCH --job-name=my_matlab_job
#SBATCH --output=my_matlab_job.out
#SBATCH --error=my_matlab_job.err
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1

module load matlab  # Load MATLAB module if necessary

# Write MATLAB script to a temporary file
echo "$MATLAB_SCRIPT" > func.m

# Loop over combinations of i and j
for i in {1..3}; do
    for j in {1..2}; do
    	for k in {1..3}; do
    	    LOG_FILE="log_${i}_${j}_${k}.txt"
        	# Submit MATLAB job to Slurm for each combination
        	srun --exclusive matlab -nodisplay -r "sequentialExptDesignBatchRunner($i, $j, $k); exit;" > "$LOG_FILE" &
        	echo "Submitted MATLAB job with ID: $JOB_ID for i=$i and j=$j and k=$k"
     done
done
