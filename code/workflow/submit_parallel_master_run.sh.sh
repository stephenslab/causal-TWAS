#!/bin/sh

#SBATCH --time=1-12:00:00
#SBATCH --array=1-20
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=40G
#SBATCH --job-name=master_run.sh
#SBATCH --output=master_run.sh_%A_%a.out
#SBATCH --error=master_run.sh_%A_%a.err
#SBATCH --partition=mstephens

head -n $SLURM_ARRAY_TASK_ID master_run.sh | tail -1 | sh
mkdir master_run.sh_JOBID_$SLURM_ARRAY_JOB_ID
