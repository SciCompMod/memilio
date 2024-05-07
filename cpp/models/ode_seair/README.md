# ODE SEAIR compartment model

This model is an extented ODE type model.
The six compartments 
- `Susceptible` ($S$), may become Exposed at any time
- `Exposed` ($E$), becomes Asymptomatic after some time
- `Asymptomatic` ($A$), becomes Infected or Recovered after some time. People which are not included in the count of confirmed cases, either because they are asymptomatic or because of insufficient testing.
- `Infected` ($I$), becomes Recovered or Dead after some time. Is assumed to be equal to confirmed cases here.
- `Recovered` ($R$)
- `Dead` ($D$)

are used to simulate the spread of the disease.

A detailed description of the model can be found in the publication
[Tsay et al. (2020), Modeling, state estimation, and optimal control for the US COVID-19 outbreak](https://doi.org/10.1038/s41598-020-67459-8).

Below is an overview of the model architecture and its compartments.
![secir_model](https://github.com/SciCompMod/memilio/assets/69154294/9d1bd9d1-9c6d-484f-aac1-6aead0b34185)

| Mathematical variable                   | C++ variable name | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\alpha_{a}$                      |  `SocialDistancing`               | Infectivity of the Asymptomatic compartment. This includes ContactPatterns, the TransmissionProbabilityOnContact and the influence of social distancing.|
| $\alpha_{i}$                      |  `Quarantined`               | Infectivity of the Infected compartment. This includes ContactPatterns, the TransmissionProbabilityOnContact and the influence of quaratine. |
| $N$                      |  `pop_total`               | Represents the total population. |
| $T_{E}$                      |  `TimeExposed`               | Time in days an individual stays in the Exposed compartment. |
| $\frac{\mu_A^R}{T_A}$                      |  `RecoveryRateFromAsymptomatic`               | The rate at which people recover from the Asymptomatic compartment. This includes the probability of not getting a positive test and the time spent in the Asymptomatic compartment before recovery.  |
| $\frac{1- \mu_A^R}{T_A}$                      |  `TestingRate`               | The rate at which people transiting from Asymptomatic to Infected compartment. This includes the probability of getting a positive test and the time spent in the Asymptomatic compartment before tested positive. |
| $\frac{\mu_I^R}{T_I}$                      |  `RecoveryRate`               | The rate at which people recover. This includes the probability of not dying and the time spent in the Infected compartment before recovery. |
| $\frac{1 - \mu_I^R}{T_I}$                     |  `DeathRate`               | The rate at which people die. This includes the probability of dying and the time spent in the Infected compartment before dying.  |
| $\frac{1}{T_R}$                      |  `TimeRecoveredInv`               | Inverse time in days an individual stays in the Recovered compartment before becoming Susceptible again. |


- An example can be found in [examples/ode_seair.cpp](../../examples/ode_seair.cpp).
- For an example of using this model to define an optimization problem, see [examples/ode_seair_optimization.cpp](../../examples/ode_seair_optimization.cpp).
