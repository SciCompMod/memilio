#!/bin/bash
# First line is the interpreter and should always be bash
# Afterwards you should at least define the job name, where the command line output should be written to

#SBATCH -J "Memilio profile" # Job name
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks-per-node=4
#SBATCH --account=training2508
#SBATCH --output=slurm/slurm-%j.out
#SBATCH --reservation=gpuhack25-day2
#SBATCH --gpus-per-task=1

#SBATCH --time=00:50:00


module load CMake GCC Eigen OpenMPI HDF5 CUDA

export NSYS_TEMP_DIR=/p/project1/training2508/MEmilio/kilian/memilio/tmp/nsys_temp_dir_$SLURM_JOB_ID
mkdir -p $NSYS_TEMP_DIR
chmod 755 $NSYS_TEMP_DIR

export TMPDIR=$NSYS_TEMP_DIR
export TMP=$NSYS_TEMP_DIR
export TEMP=$NSYS_TEMP_DIR

echo "Temp Dir: $NSYS_TEMP_DIR" 


# srun --cpu-bind=core nsys profile -o $NSYS_TEMP_DIR/profile_%h_%p ./bin/ode_secir_graph_example

nsys profile -o $NSYS_TEMP_DIR/profile_new_${SLURM_JOB_ID} --trace nvtx,cuda,openacc $PROJECT_training2508/MEmilio/kilian/cpp/build/bin/ode_secir_graph_example 
#$PROJECT_training2508/MEmilio/kilian/cpp/build/bin/ode_secir_graph_example
