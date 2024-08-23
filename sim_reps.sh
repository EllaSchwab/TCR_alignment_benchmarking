#!/bin/bash
#SBATCH --job-name=immuneSIM
#SBATCH --output=immuneSIM.out
#SBATCH --error=immuneSIM.err
#SBATCH --time=05:00:00  # Adjust the time limit as needed
#SBATCH --mem=8G
#SBATCH --partition=main  # Adjust the partition/queue as needed


module purge
module load gcc/11.3.0
module load openblas/0.3.20
module load r/4.4.0



# Run the R script
Rscript GenerateRepertoires.R
