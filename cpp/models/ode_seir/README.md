# ODE SEIR compartment model

This model is a very simple ODE model with only four compartments and few parameters, mostly for demonstration of the MEmilio framework:
- Susceptible, may become exposed at any time
- Exposed, becomes infected after some time
- Infected, will recover after some time
- Recovered

Below is an overview of the model architecture and its compartments.

![SEIR_model](https://github.com/DLR-SC/memilio/assets/69154294/11ee01be-40dc-40d0-9157-1a4bec775b02)
| Parameter                   | Implementation | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\phi$                      |  `ContactPatterns`               | Daily contact rate |
| $\rho$                      |  `TransmissionProbabilityOnContact`               | Transmission risk for people located in one of the susceptible compartments |
| $N$                         | `populations.get_total()`   | Total population. |
| $T_{E}$                    |  `TimeExposed`               | Time in days an individual stays in the Exposed compartment. |
| $T_{I}$                    |  `TimeInfected`               | Time in days an individual stays in the Infected compartment. |


An example can be found in [examples/ode_seir.cpp](../../examples/ode_seir.cpp)