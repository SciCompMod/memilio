#!/bin/bash
#SBATCH --job-name=lct-performance
#SBATCH --output=lct-%A.out
#SBATCH --error=lct-%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --exclusive
#SBATCH --exclude="be-cpu05, be-gpu01"
#SBATCH --time=5-0:00:00

## This script can be used to monitor runtimes for the lct model using the file lct_runtime.cpp.
## The command "sbatch get_runtimes_lct.sh ./bin/lct_runtime" should be run in this folder 
## to have valid folder instructions.

warm_up_runs=10
runs=100
echo Running $1 on node $SLURM_JOB_NODELIST with $warm_up_runs warm up runs and $runs runs.

cd ../../
if [ ! -d "build/" ]; then
    mkdir "build/"
fi
cd build/
cmake -DCMAKE_BUILD_TYPE="Release" -DMEMILIO_ENABLE_OPENMP=ON ..

for i in {1..200}
do
    cmake -DNUM_SUBCOMPARTMENTS=$i -DCMAKE_BUILD_TYPE="Release" -DMEMILIO_ENABLE_OPENMP=ON .
    cmake --build . --target lct_runtime
    srun --cpus-per-task=1 --cpu-bind=cores ./$1 $warm_up_runs $runs
done
