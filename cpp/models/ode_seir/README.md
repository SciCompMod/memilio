# ODE SEIR compartment model

This model is a very simple ODE model with only four compartments and few parameters, mostly for demonstration of the MEmilio framework:
- Susceptible, may become exposed at any time
- Exposed, becomes infected after some time
- Infected, will recover after some time
- Recovered

Below is an overview of the model architecture and its compartments.

![SEIR_model](https://github.com/DLR-SC/memilio/assets/69154294/80a36be5-57d9-4012-9b5f-25eb08ec8837)
| Mathematical variable                   | C++ variable name | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\phi$                      |  `ContactPatterns`               | Daily contact rate / Number of daily contacts. |
| $\rho$                      |  `TransmissionProbabilityOnContact`               | Transmission risk for people located in the Susceptible compartment. |
| $N$                         | `populations.get_total()`   | Total population. |
| $T_{E}$                    |  `TimeExposed`               | Time in days an individual stays in the Exposed compartment. |
| $T_{I}$                    |  `TimeInfected`               | Time in days an individual stays in the Infected compartment. |


An example can be found in [examples/ode_seir.cpp](../../examples/ode_seir.cpp)