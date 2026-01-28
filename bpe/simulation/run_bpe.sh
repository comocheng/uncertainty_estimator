#!/bin/bash
#SBATCH --job-name=mk3
#SBATCH --error=error.log
#SBATCH --nodes=1
#SBATCH --partition=short,west
#SBATCH --mem=80Gb
#SBATCH --time=4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=24


# python run_bpe.py
mpiexec -n $SLURM_NTASKS -v python run_bpe.py
