#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --exclusive
#SBATCH -t 5-0:00:00
#SBATCH --output=steps_equationbased-%A.out
#SBATCH --error=steps_equationbased-%A.err
#SBATCH --exclude="be-cpu05, be-gpu01"
#SBATCH --job-name=steps_mobilitymodels

echo Running $1 on node $SLURM_JOB_NODELIST.
cd ../cpp/build
cmake -DCMAKE_BUILD_TYPE="Release" -DMEMILIO_ENABLE_OPENMP=ON ..
cmake --build . --target $1

for i in $(seq 2 12) 
do
    srun --cpu-bind=core --cpus-per-task=1 ./bin/$1 1e-$i
done