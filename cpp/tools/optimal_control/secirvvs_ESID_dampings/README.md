# ESID Damping Optimization

This tool solves an optimal control problem for the SECIRVVS / ESID epidemic model with time-dependent damping interventions.

It reads intervention definitions and constraint bounds, computes an optimal intervention schedule, and writes the result as CSV files.

# Usage

`./optimal_control_ESID <data_dir> <settings_dir>`

The settings directory is already provided in the tools folder `/optimal_control_settings/` and contains: 
- intervention_list.json
- constraints.json

The data directory is used to create the SECIRVVS model.

# Tool layout

The program entry point is the script
`secirvvs_ESID_dampings/optimal_control_ESID.cpp`.
It loads the inputs, runs the optimization, writes results.

We provide a brief overview of the toolbox structure.

### Optimal control configuration: `optimal_control_settings/`
- `intervention_list.json`
List of interventions to be used in the optimization.
- `constraints.json`
Path and terminal constraints on epidemic states.

### Create constraints: `create_constraints_ESID.h`
- Load the path- and terminal constraint from file (constraints.json)

### Define control parameters: `create_control_parameters.h`
- Defines control parameters (damping, costs, effectiveness)
- Loads interventions from file (intervention_list.json)

### Check feasability `check_feasability.h`
- Run a simulation with the most restrictive parameters
- If the constraints are not fulfilled -> infeasible problem

### Simulation model used in optimization: `optimization_model/`
- Builds the SECIRVVS model
- Uses a synthetic setup if no real data is loaded

### Optimization configuration: `optimization_settings/`
- Stores settings and configurations of everything

### Solver interface: `ipopt_solver/`
- Defines objective and constraint evaluations
- Connects simulations to Ipopt
- Saves final results

### Small utilities: `helpers/`
- Time grid creation
- Integrator selection
- Numerical helpers

### Define application of controls: `control_parameters/`
- Defines the class to store the control parameters

### Define handling of constraints: `constraints/`
- Maps constraint names to model states
- Evaluates path and terminal constraints during simulation

-----

# Outputs

Written to the working directory:

- `control_parameters.csv`
Time series of optimized intervention strengths.
- `population_time_series.csv`
Epidemic trajectory under optimal controls.
- `constraint_violation.log`
Only written if the feasibility check fails.

# Plotting

- `python3 plot_controls.py` 
- `python3 plot_secirvvs.py`

# Notes

- Feel free to change the model to an actual one instead of an artificial example
- You must use a controlled integrator (ControlledCashKarp54 or ControlledFehlberg78) with a integrating resolution
of at least >=10 or 20 to be safe. The correct gradient evaluation requires higher precision than standard values.
- Reverse mode requires multiple GB of memory which can result in segmentation faults.
- We use forward mode since it is faster when lots of cores are available.
- Activation is linear since we allow for continous controls. We don't want to push controls to integer values (0 or 1).
- Maximilian Betz created a reduced version with a simpler file strucutre for easier ESID integration
