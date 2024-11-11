#!/bin/bash
#SBATCH --job-name=lct-cache
#SBATCH --output=lct-%A.out
#SBATCH --error=lct-%A.err
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --exclusive
#SBATCH --exclude="be-cpu05, be-gpu01"
#SBATCH -t 5-0:00:00

warm_up_runs=0
runs=10
echo Running $1 on node $SLURM_JOB_NODELIST with $warm_up_runs warm up runs and $runs runs.
cd ../build
for i in {1..50}
do
    cmake -DNUM_SUBCOMPARTMENTS=$i .
    cmake --build . --target lct_timing
    srun --cpu-bind=core valgrind --tool=callgrind --simulate-cache=yes ./$1  $warm_up_runs $runs
done


