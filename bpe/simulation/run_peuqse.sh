#!/bin/bash
#SBATCH --job-name=PEUQSE
#SBATCH --error=error.log
#SBATCH --nodes=1
#SBATCH --partition=short,west
#SBATCH --mem=20Gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=24


python run_peuqse.py

