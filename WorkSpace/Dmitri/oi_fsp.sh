#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=3:00:00
#SBATCH --partition=csu
#SBATCH --qos=normal
#SBATCH --ntasks=32
#SBATCH --job-name=OI_FSP
#SBATCH --output=OI_FSP.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=svetlov@colostate.edu

module purge

module load matlab/R2023b

matlab -nodesktop -nodisplay -r 'clear; OI_FSP;'