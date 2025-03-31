#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --exclusive
#SBATCH -t 5-0:00:00
#SBATCH --output=likwid_graphbased-%A.out
#SBATCH --error=likwid_graphbased-%A.err
#SBATCH --exclude="be-cpu05, be-gpu01"
#SBATCH --job-name=likwid_mobilitymodels

warm_up_runs=0
runs=1
regions=80

module purge
module load PrgEnv/gcc13-openmpi

echo Running $1 on node $SLURM_JOB_NODELIST with $warm_up_runs warm up runs and $runs runs.
cd ../cpp/build
rm -rf CMakeCache.txt CMakeFiles/
CXX="g++ -fverbose-asm" cmake -DCMAKE_BUILD_TYPE="Release" -DMEMILIO_ENABLE_OPENMP=ON -DMEMILIO_ENABLE_WARNINGS_AS_ERRORS=OFF ..
# cmake --build . --target graph_timing

# perf record ./bin/graph_timing $warm_up_runs $runs $regions
# perf report 
srun --cpu-bind=cores --cpus-per-task=1 --cpu-freq=2200000-2200000 likwid-perfctr -C 0 -g MEM_DP -m ./bin/graph_timing $warm_up_runs $runs $regions
