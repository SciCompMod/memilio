# SDE-based SEIR-type model with two variants

This module models and simulates the epidemic using an SDE-based SEIR-type model approach with two variants. After a recovered infection from the first variant (considered the wild type) you can still get infected with the second variant. Infection with the second variant immunizes you against the first variant. Unlike the agent based model that uses particular agents, this model simulates the spread of a communicable disease in a population with subpopulations being in different compartments. The compartments are as follows: `Susceptible`,  `ExposedV1`, `ExposedV2`, `ExposedV1V2`,  `InfectedV1`, `InfectedV2`, `InfectedV1V2`, `RecoveredV1`, `RecoveredV2` and `RecoveredV1V2`. The compartments are addended by the relevant variants. The addendum `V1` means that the individual is infected with the first variant, the addendum `V2`means that the individual is infected with the second variant with no prior infections and the addendum `V1V2` means that the individual is infected with the second variant after a successful recovery from the first variant. Only individuals in the infected compartments are infectious. 




# SDE SEIR two variant compartment model

This model is a stochastic multivariant SDE SEIR model with  10 compartments and few parameters, addressing reinfection with a second variant. 
- Susceptible, may become exposed at any time.
- Exposed, may become infectious at any time
- Infected, will be recovered after some time.
- Recovered, recovered from previous infection, 

We assume simulations over short periods of time, so that the population size can be considered constant and birth as well as (natural) mortality rates can be ignored. 

Important note on the solution: The numerical integration method is a Euler-Maruyama which reuses the EulerIntegratorCore of MEmilio accordingly. The (re)use of other implemented integrators for ODEs is neither considered nor suggested at the moment.



