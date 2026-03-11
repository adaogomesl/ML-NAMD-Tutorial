#!/bin/sh
## script for PyRAI2MD
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-00:00:00
#SBATCH --job-name=dbh-1
#SBATCH --partition=lopez,short,long
#SBATCH --mem=110mb
#SBATCH --output=%j.o.slurm
#SBATCH --error=%j.e.slurm

export INPUT=input
export WORKDIR=/projects/lopez/share_from_Leticia/ml-nac/NEQUIP/on-paper-models-train-before-hpo/dbh-4traj-and-90b-60a/4traj-and-90b-60a-pc/2026-01-23-16-55-49-946501/500-namd/dbh-1

cd $WORKDIR
pyrai2md $INPUT

