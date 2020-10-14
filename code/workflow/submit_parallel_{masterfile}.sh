#!/bin/sh

#SBATCH --time=1-12:00:00
#SBATCH --array=1-wc:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --job-name={masterfile}
#SBATCH --output={masterfile}_%A_%a.out
#SBATCH --error={masterfile}_%A_%a.err
#SBATCH --partition=mstephens

head -n $SLURM_ARRAY_TASK_ID {masterfile} | tail -1 | sh
mkdir {masterfile}_JOBID_$SLURM_ARRAY_JOB_ID
