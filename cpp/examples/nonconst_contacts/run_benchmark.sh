#!/bin/bash
#SBATCH --job-name=smoother-runtime
#SBATCH --output=smoother-runtime%A.out
#SBATCH --error=smoother-runtime%A.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 90
#SBATCH --exclusive
#SBATCH --exclude="be-cpu05, be-gpu01"

# Run this script from build folder with downloaded data in repository
num_runs=100
num_warm_up_runs=20

# Load module
module purge
module load PrgEnv/gcc13-openmpi

srun --cpus-per-task=1 --cpu-bind=core ./bin/runtime_comparison_smoother -NumberRuns $num_runs -NumberWarmupRuns $num_warm_up_runs
