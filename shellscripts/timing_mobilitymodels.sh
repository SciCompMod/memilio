#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --exclusive
#SBATCH -t 5-0:00:00
#SBATCH --output=timing_equationbased_noage_euler-%A.out
#SBATCH --error=timing_equationbased_noage_euler-%A.err
#SBATCH --exclude="be-cpu05, be-gpu01"
#SBATCH --job-name=timing_mobilitymodels

warm_up_runs=10
runs=100

echo Running $1 on node $SLURM_JOB_NODELIST with $warm_up_runs warm up runs and $runs runs.
cd ../cpp/build
rm -rf CMakeCache.txt CMakeFiles/
cmake -DCMAKE_BUILD_TYPE="Release" -DMEMILIO_ENABLE_OPENMP=ON ..
cmake --build . --target $1
for i in $(seq 1 1 4) 
do
    srun --cpu-bind=core --cpus-per-task=1 ./bin/$1 $warm_up_runs $runs $i
done
for i in $(seq 300 5 350) 
do
    srun --cpu-bind=core --cpus-per-task=1 ./bin/$1 $warm_up_runs $runs $i
done
