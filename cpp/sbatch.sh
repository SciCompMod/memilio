#!/bin/bash -x
#SBATCH --account=training2508
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --output=out/slurm-%j-acc.out
#SBATCH --error=out/slurm-%j-acc-err
#SBATCH --time=00:15:00
#SBATCH --partition=all

lscpu

module load CMake GCC Eigen OpenMPI HDF5 CUDA

srun --gpufreq=1410 ./build/bin/ide_secir_hackathon_example
srun --gpufreq=1410 ./build/bin/ide_secir_hackathon_acc_example

# nsys profile ./build/bin/ide_secir_hackathon_acc_example