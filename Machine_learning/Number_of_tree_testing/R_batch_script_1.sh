#!/bin/bash

#SBATCH -p defq
#SBATCH --job-name=R_ML_1-200
#SBATCH --output=Rscript_1-200-%J.txt

module load conda/3.6
source activate R_env 

Rscript ML_number_trees_1-200.R

