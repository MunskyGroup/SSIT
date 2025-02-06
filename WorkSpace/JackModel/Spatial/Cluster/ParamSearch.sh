#!/bin/bash
#SBATCH --job-name=matlab_parallel_job   # Name of the job
#SBATCH --output=matlab_parallel_output_%j.txt   # Output file (%j will be replaced with the job ID)
#SBATCH --error=matlab_parallel_error_%j.txt    # Error file (%j will be replaced with the job ID)
#SBATCH --ntasks=1                    # Number of tasks (usually 1 for MATLAB)
#SBATCH --cpus-per-task=4             # Number of CPUs per task (adjust based on parallelization needs)
#SBATCH --mem=16G                     # Amount of memory requested (adjust as needed)
#SBATCH --time=02:00:00               # Maximum runtime in the format HH:MM:SS
#SBATCH --partition=compute           # Partition or queue to use (adjust to your cluster's configuration)

# Load MATLAB module (adjust to your cluster's module system)
module load matlab

# Run the MATLAB script
matlab -nodisplay -r "try, ParamSearch, catch ME, disp(ME.message), exit(1), end, exit(0)"