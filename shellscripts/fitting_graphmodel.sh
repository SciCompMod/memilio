#!/bin/bash 
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 5-0:00:00
#SBATCH --output=shellscripts/create_testdata_fitting-%A.out
#SBATCH --error=shellscripts/create_testdata_fitting-%A.err
#SBATCH --exclude="be-cpu05, be-gpu01"
#SBATCH --job-name=create_testdata_fitting

module load PrgEnv/gcc13-openmpi-python
source venv/bin/activate
srun --cpu-bind=core python pycode/examples/simulation/graph_germany_nuts3.py