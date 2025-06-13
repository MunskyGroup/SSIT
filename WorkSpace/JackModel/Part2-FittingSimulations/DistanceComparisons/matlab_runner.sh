#!/bin/bash
#SBATCH --job-name=matlab_job   # Name of the job

# Load MATLAB module (adjust to your cluster's module system)
module load matlab-2022b
module load gnu9/9.4.0

# Start timing
start_time=$(date +%s)

fileLocation=$1
args=$2

mkdir -p "${PWD}/output"
output_names="${PWD}/output/output_${SLURM_JOB_ID}_$(basename ${fileLocation})"

if [ -z "$args" ]; then
    matlab -nodisplay -r "try, run('${fileLocation}'), catch ME, disp(ME.message), end, exit;" \
    >> "$output_names" 2>&1 &
else
    matlab -nodisplay -r "try, ${fileLocation}(${args}), catch ME, disp(ME.message), end, exit;" \
    >> "$output_names" 2>&1 &
fi

wait

end_time=$(date +%s)
total_time=$(( (end_time - start_time) / 60 ))
echo "Total time to complete the job: $total_time minutes"
exit 0
