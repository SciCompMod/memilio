# ODE SEAIR compartment model

This model is an extented ODE type model.
The six compartments 
- `Susceptible` ($S$), may become Exposed at any time
- `Exposed` ($E$), becomes Asymptomatic after some time
- `Asymptomatic` ($A$), becomes Infected or Recovered after some time
- `Infected` ($I$), becomes Recovered or Dead after some time
- `Recovered` ($R$)
- `Dead` ($D$)

are used to simulate the spread of the disease.

A detailed description of the model can be found in the publication
[Tsay et al. (2020), Modeling, state estimation, and optimal control for the US COVID-19 outbreak](https://doi.org/10.1038/s41598-020-67459-8). 

| Mathematical variable                   | C++ variable name | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\alpha_{a}$                      |  ``               | . |
| $\alpha_{i}$                      |  ``               | . |


- An example can be found in [examples/ode_seair.cpp](../../examples/ode_seair.cpp).
- For an example of using this model to define an optimization problem, see [examples/ode_seair_optimization.cpp](../../examples/ode_seair_optimization.cpp).
