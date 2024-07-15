# SDE-based SEIR-type model with two variants

This module models and simulates the epidemic using an SDE-based SEIR-type model approach with two variants. After a recovered infection from the first variant (considered the wild type) you can still get infected with the second variant. Infection with the second variant immunizes you against the first variant. Unlike the agent based model that uses particular agents, this model simulates the spread of a communicable disease in a population with subpopulations being in different compartments. The model also introduces stoachsticity compared to an ODE-model. The compartments are as follows: `Susceptible`,  `ExposedV1`, `ExposedV2`, `ExposedV1V2`,  `InfectedV1`, `InfectedV2`, `InfectedV1V2`, `RecoveredV1`, `RecoveredV2` and `RecoveredV1V2`. The compartments are addended by the relevant variants. The addendum `V1` means that the individual is infected with the first variant, the addendum `V2`means that the individual is infected with the second variant with no prior infections and the addendum `V1V2` means that the individual is infected with the second variant after a successful recovery from the first variant. Only individuals in the infected compartments are infectious. 

Below is an overview of the model architecture.



## Simulation

The simulation runs in discrete time steps using an euler-maruyama integration scheme. The Simulation class handles the parameters and the numerical integrator. It also stores the result. 
| Mathematical variable                   | C++ variable name | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\phi$                      |  `ContactPatterns`               | Daily contact rate / Number of daily contacts. |
| $\rho_1$                      |  `TransmissionProbabilityOnContactV1`               | Transmission risk for people located in the Susceptible compartment (susceptible to infection with variant 1). |
| $\rho_2$                      |  `TransmissionProbabilityOnContactV2`               | Transmission risk for people located in the Susceptible compartment or in the RecoveredV1 compartment (susceptivle to infection with variant 2). |
| $N$                         | `populations.get_total()`   | Total population. |
| $T_{E1}$                    |  `TimeExposed1`               | Average time in days an individual stays in the ExposedV1 compartment. |
| $T_{E2}$                    |  `TimeExposed2`               | Average time in days an individual stays in the ExposedV2 or in the ExposedV1V2 compartment. |
| $T_{I1}$                    |  `TimeInfected1`               | Average time in days an individual stays in the InfectedV1 compartment. |
| $T_{I1}$                    |  `TimeInfected2`               | Average time in days an individual stays in the InfectedV1 or in the InfectedV2 compartment. |

An example can be found in [examples/sde_seirvv.cpp](../../examples/sde_seirvv.cpp)