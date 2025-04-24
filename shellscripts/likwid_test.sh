#!/bin/sh
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH --exclusive
#SBATCH -t 5-0:00:00
#SBATCH --output=likwid_test-%A.out
#SBATCH --error=likwid_test-%A.err
#SBATCH --exclude="be-cpu05, be-gpu01"
#SBATCH --job-name=likwid_test


cd ../cpp/examples
g++ -O3 -o likwid_test likwid_test.cpp

srun --cpu-bind=cores --cpus-per-task=1 --cpu-freq=2200000-2200000 likwid-perfctr -C S0:1 -g MEM_DP ./likwid_test