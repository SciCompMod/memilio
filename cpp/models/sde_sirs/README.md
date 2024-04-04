
# SDE SIRS compartment model

This model is a simple stochastic SDE SIRS model with only three compartments and few parameters, addressing waning immunity and allowing reinfection after some time of immunity.
- Susceptible, may become infected at any time.
- Infected, will be recovered after some time.
- Recovered, recovered from previous infection, temporarily immune.

We assume simulations over short periods of time, so that the population size can be considered constant and birth as well as (natural) mortality rates can be ignored. 

Important note on the solution: The numerical integration method is a Euler-Maruyama which reuses the EulerIntegratorCore of MEmilio accordingly. The (re)use of other implemented integrators for ODEs is neither considered nor suggested at the moment.
Below is an overview of the model architecture and its compartments.

![SIR_model](https://github.com/SciCompMod/memilio/assets/69154294/bdde5ca3-cd7d-4695-8704-03a712c92ff7)
| Mathematical variable                   | C++ variable name | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\phi$                      |  `ContactPatterns`               | Daily contact rate / Number of daily contacts. |
| $\rho$                      |  `TransmissionProbabilityOnContact`               | Transmission risk for people located in the Susceptible compartment. |
| $N$                         | `populations.get_total()`   | Total population. |
| $T_{I}$                    |  `TimeInfected`               | Time in days an individual stays in the Infected compartment. |
| $T_{R}$                    |  `TimeImmune`               | Time in days an individual stays in the Recovered compartment before becoming Susceptible again. |

An example can be found in [examples/ode_sir.cpp](../../examples/ode_sir.cpp)
