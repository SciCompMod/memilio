# Scaling examples

Build with the following command from cpp folder:

```bash
mkdir build && cd build
cmake ..
cmake --build . --target ode_runtime lct_runtime ide_runtime ode_ensemble_runs lct_ensemble_runs ide_ensemble_runs
```

To run the code use the .sh files "get_runtime.sh" and "get_ensemble_runs.sh" in the three folders "simulation_paper_*".
The folder "simulation_paper" is not part of the scaling examples.
Start the batch files with the current directory beeing the build folder (this is the case after following the build commands).

You can find an [runtime estimate](../plots_paper/calculate_runtime.py) for all 6 examples.

## Scaling of graph

Look at runtime increase with an increasing number of nodes in a graph.
ODE and LCT have 10% mobile population evenly distributed for each node to its 100 neighboors.
IDE does not use mobility.

- ode_runtime
- lct_runtime
- ide_runtime

Results are in the folder "results_runtime" in each of the model paper example subfolder.

## Parallel scaling

Look at runtime of ParameterStudy with an increasing number of cpu nodes to distribute on.
Number of runs and model is the same, while some parameters include uncertainty. 
Models use rki data to be close to real world use cases.

ODE uses a graph with 400 nodes for german county. Due to high runtime per graph, the number of runs is only 1280 (10 runs for highest cpu number that is multiple of 2).
IDE and LCT do implement a single model for germany. Thus number of runs is 16384 (divisible by 128)

- ode_ensemble_runs
- lct_ensemble_runs
- ide_ensemble_runs

Results are on the terminal (due to the batch script inside a file that saves the terminal output).

## Plots

Can be found in "plots_paper"
