#!/bin/bash
#
#SBATCH --job-name=compgen
#SBATCH --output=/data/douglaslab/douglste/script_logs/slurm-%A.out
#SBATCH --account=douglaslab
#SBATCH --partition=douglaslab,node
#
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=16gb
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=douglste@lafayette.edu

srun hostname

source ~/.bashrc

source activate molusc

cd ~/projects/MOLUSC/

srun python scripts/generate_companions.py
