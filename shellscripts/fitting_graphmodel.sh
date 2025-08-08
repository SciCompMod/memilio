#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 5-0:00:00
#SBATCH --output=shellscripts/train_diffusion_model-%A.out
#SBATCH --error=shellscripts/train_diffusion_model-%A.err
#SBATCH --job-name=train_diffusion_model
#SBATCH --partition=gpu
#SBATCH --gpus=1

module load PrgEnv/gcc13-openmpi-python
module load cuda/12.9.0-none-none-6wnenm2
source venv/bin/activate
srun --cpu-bind=core python pycode/examples/simulation/graph_germany_nuts3.py   