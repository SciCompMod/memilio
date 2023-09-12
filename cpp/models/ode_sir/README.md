
# ODE SIR compartment model

This model is a very simple ODE model with only three compartments and few parameters, mostly for demonstration of the MEmilio framework:
- Susceptible, may become infected at any time
- Infected, will be removed after some time
- Removed, removed from infectious process (dead or recovered).
We assume simulations over short periods of time, so that the population size can be considered constant and birth as well as (natural) mortality rates can be ignored. 

An example can be found in [examples/ode_sir.cpp](../../examples/ode_sir.cpp)
