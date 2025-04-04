#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --exclusive
#SBATCH -t 5-0:00:00
#SBATCH --output=likwid_equationbased-%A.out
#SBATCH --error=likwid_equationbased-%A.err
#SBATCH --exclude="be-cpu05, be-gpu01"
#SBATCH --job-name=likwid_mobilitymodels

warm_up_runs=0
runs=50
regions=400

module purge
module load PrgEnv/gcc12-openmpi

echo Running $1 on node $SLURM_JOB_NODELIST with $warm_up_runs warm up runs and $runs runs.
cd ../cpp/build
rm -rf CMakeCache.txt CMakeFiles/
cmake -DCMAKE_BUILD_TYPE="Release" -DMEMILIO_ENABLE_OPENMP=ON -DMEMILIO_ENABLE_WARNINGS_AS_ERRORS=OFF ..
cmake --build . --target ode_metapop_timing

perf record -C 0 likwid-pin ./bin/ode_metapop_timing $warm_up_runs $runs $regions
perf report 
# srun --cpu-bind=cores --cpus-per-task=1 --cpu-freq=2200000-2200000 valgrind --tool=massif --detailed-freq=2 ./bin/ode_metapop_timing $warm_up_runs $runs $regions
# srun --cpu-bind=cores --cpus-per-task=1 --cpu-freq=2200000-2200000 likwid-perfctr -C 0 -g MEM_DP -m ./bin/ode_metapop_timing $warm_up_runs $runs $regions
