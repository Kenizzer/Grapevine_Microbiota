#!/bin/bash

#SBATCH -p defq
#SBATCH --job-name=R_ML_201-400
#SBATCH --output=Rscript_201-400-%J.txt

module load conda/3.6
source activate R_env 

Rscript ML_number_trees_201-400.R



