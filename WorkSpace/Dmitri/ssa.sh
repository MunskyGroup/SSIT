#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --partition=csu
#SBATCH --qos=normal
#SBATCH --ntasks=16
#SBATCH --job-name=NFDP_SSA
#SBATCH --output=NFDP_SSA.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=svetlov@colostate.edu

module purge

module load matlab/R2023b

matlab -nodesktop -nodisplay -r 'clear; NFDP_SSA;'