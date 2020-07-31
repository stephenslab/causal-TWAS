#!/bin/sh

#SBATCH --time=1-12:00:00
#SBATCH --array=1-wc:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --job-name=master_run1-x.sh
#SBATCH --output=master_run1-x.sh_%A_%a.out
#SBATCH --error=master_run1-x.sh_%A_%a.err
#SBATCH --partition=mstephens

head -n $SLURM_ARRAY_TASK_ID master_run1-x.sh | tail -1 | sh
mkdir master_run1-x.sh_JOBID_$SLURM_ARRAY_JOB_ID
