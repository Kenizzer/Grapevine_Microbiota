#!/bin/bash

#SBATCH -p defq
#SBATCH --job-name=Jup-notebook-job
#SBATCH --output=jupyter-log-%J.txt

module load conda/3.6
source activate R_env 


## get tunneling info
XDG_RUNTIME_DIR=""
ipnport=8999
ipnip=$(hostname -i)


## print tunneling instructions to jupyter-log-{jobid}.txt
echo -e "
    Copy/Paste this in your local terminal to ssh tunnel with remote
    -----------------------------------------------------------------
    ssh -N -L $ipnport:$ipnip:$ipnport user@host
    -----------------------------------------------------------------

    Then open a browser on your local machine to the following address
    ------------------------------------------------------------------
    localhost:$ipnport  (prefix w/ https:// if using password)
    ------------------------------------------------------------------
    "

## start an ipcluster instance and launch jupyter server
jupyter-notebook --no-browser --port=$ipnport --ip=$ipnip

source deactivate
