#!/bin/bash

#SBATCH -p defq
#SBATCH --job-name=R_ML_401-501
#SBATCH --output=Rscript_401-501-%J.txt

module load conda/3.6
source activate R_env 

Rscript ML_number_trees_401-501.R



