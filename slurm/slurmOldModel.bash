#!/bin/bash
#SBATCH -J old-model                # job name
#SBATCH -o old-model.o%j            # output and error file name (%j expands to SLURM jobID)
#SBATCH -N 1                        # number of nodes requested
#SBATCH -n 128                      # total number of tasks to run in parallel
#SBATCH -p normal                   # queue (partition) 
#SBATCH -t 48:00:00                 # run time (hh:mm:ss) 
#SBATCH -A MCB24023                 # Allocation name to charge job against

module load launcher
source /work/10081/kellyktvt/ls6/pinetree/.venv-trnas/bin/activate

export LAUNCHER_WORKDIR=/work/10081/kellyktvt/ls6/phage_simulation
export LAUNCHER_JOB_FILE=/work/10081/kellyktvt/ls6/phage_simulation/src/python/models/jobs_phage_model_fixed.txt

${LAUNCHER_DIR}/paramrun