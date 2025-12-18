#!/bin/bash
#SBATCH --job-name=ode-runtime
#SBATCH --output=ode-runtime%A.out
#SBATCH --error=ode-runtime%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=17
#SBATCH --cpus-per-task=1
#SBATCH -t 90
#SBATCH --exclusive
#SBATCH --exclude="be-cpu05, be-gpu01"

# Run this script from build folder with downloaded data in repository
num_runs=100
num_warm_up_runs=10

# Load module
module purge
module load PrgEnv/gcc13-openmpi

for i in {200..1000..50}
do
    echo Run with $i regions.
    srun --cpus-per-task=1 --cpu-bind=core ./bin/ode_runtime -NumberRuns $num_runs -NumberWarmupRuns $num_warm_up_runs -NumberRegions $i
done