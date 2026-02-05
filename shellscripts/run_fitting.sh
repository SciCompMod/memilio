#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 5-0:00:00
#SBATCH --output=shellscripts/run_fitting-%A.out
#SBATCH --error=shellscripts/run_fitting-%A.err
#SBATCH --job-name=run_fitting

module load PrgEnv/gcc13-openmpi-python
module load cuda/12.9.0-none-none-6wnenm2

file=graph_spain_nuts3.py
source venv/bin/activate
echo "Running script: $file"
srun --cpu-bind=core python pycode/examples/simulation/$file