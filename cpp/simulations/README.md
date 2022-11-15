# Simulations #

This folder provides different simulations of spatially resolved graph-ODE (or metapopulation) models
using one ODE model for each county and realizing inter-county mobility via a graph approach.

- 2020_npis_wildtype: Focus on a SECIR model using parameters for Sars-CoV-2 wild type variant and
implementing static nonpharmaceutical interventions (NPIs) as well as dynamic NPIs. Dynamic NPIs
get into play once predefined incidence thresholds (50 and 200) are exceeded.

- 2021_vaccination_delta: Extending the model of 2020_npis_wildtype by vaccination and reinfection and
considering the effect of vaccination in combination with the lifting of NPIs during the arrival of Delta.

In the following we will show how to run the simulations.

## Run instruction

Here we will use the  ``2021_vaccination_sarscov2_delta_germany`` simulation as example to show how to execute them. The application for the other simulation is analogous.

In general, the simulations requieres a set of data provided in the [MEmilio Epidata Package](#https://github.com/DLR-SC/memilio/tree/main/pycode/memilio-epidata).
The way we can run the files there is well explained. Once we have activated the virtual environment and installed the Epidata package, we start generating the data.
For all simulations, we use a seven day moving average. So before running all python scripts for data generation, we need to set the moving average.
The simplest option are either to set the default value in getDataIntoPandasDataFrame.py to 7 or to ``add arg_dict['moving_average'] = 7`` in all mains.
Once we done that, we should run all the python files beginning with ``get..Data``. When this is done, a lot of data has been stored in our data folder.

Next, we start the simulation. There are several different possibilities to run this. 
For example, we can read in a graph or create a new one and also specify different run options.
More details about the run options can be found in the simulation file. Here we will only show the basic application.
For the simplest functionality, we only need the path to the data folder, one save folder and one results folder.


The simulation can then be started with these arguments:

```C++
/localdata1/code/memilio$ cd cd build/simulations/
/localdata1/code/memilio/build/simulations$ ./2021_vaccination_delta "/localdata1/code/memilio/data" "save_dir" "results_dir" 
```
The number of runs and simulation days can be modified in the run function.

When we have successfully calculated the simulation, the results can still be visualized. For this we refer here to the [tools](#https://github.com/DLR-SC/memilio/tree/main/tools) folder

