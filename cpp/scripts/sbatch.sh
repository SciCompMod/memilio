#!/bin/bash
#SBATCH --account=training2508
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cores=64
#SBATCH --threads-per-core=1
#SBATCH --gres=gpu:4
#SBATCH --gpus-per-task=1
#SBATCH --output=slurm-%j-acc.out
#SBATCH --error=slurm-%j-acc-err
#SBATCH --time=00:15:00
#SBATCH --partition=all
#SBATCH --reservation gpuhack25-day4

module load NVHPC CMake CUDA NVHPC Eigen HDF5 OpenMPI

export OMP_NUM_THREADS=64

# srun --cpus-per-task=$OMP_NUM_THREADS nsys profile --trace cuda,openacc,nvtx ./../build/bin/ide_secir_example


srun --cpus-per-task=64 ./../build/bin/ide_secir_example

# srun nsys profile --trace cuda,openacc,nvtx ./../build/bin/ide_secir_example



# srun: error: Invalid numeric value "" for --cpus-per-task.
