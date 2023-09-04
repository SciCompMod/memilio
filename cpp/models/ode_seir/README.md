# ODE SEIR compartment model

This model is a very simple ODE model with only four compartments and few parameters, mostly for demonstration of the MEmilio framework:
- Susceptible, may become exposed at any time
- Exposed, becomes infected after some time
- Infected, will recover after some time
- Recovered

An example can be found in [examples/ode_seir.cpp](../../examples/ode_seir.cpp)
