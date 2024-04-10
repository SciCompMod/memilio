# IDE SECIR model

This model is based on Integro-differential equations.
The eight compartments 
- `Susceptible` ($S$), may become exposed at any time
- `Exposed` ($E$), becomes infected after some time
- `InfectedNoSymptoms` ($I_{NS}$), becomes InfectedSymptoms or Recovered after some time
- `InfectedSymptoms` ($I_{Sy}$), becomes InfectedSevere or Recovered after some time
- `InfectedSevere` ($I_{Sev}$), becomes InfectedCritical or Recovered after some time
- `InfectedCritical` ($I_{Cr}$), becomes Recovered or Dead after some time
- `Recovered` ($R$)
- `Dead` ($D$)

are used to simulate the spread of the disease.

Below is an overview of the model architecture and its compartments.

![tikzIDESECIR](https://github.com/SciCompMod/memilio/assets/70579874/3500421a-035c-4ce1-ae95-a54d8097be82)

The model parameters used are the following:

| Mathematical variable                   | C++ variable name | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\phi$                      |  `ContactPatterns`               | Average number of contacts of a person per day. |
| $k$                      |  `Seasonality`               |  The influence of the seasons is taken into account with the seasonality parameter. |
| $\rho$                      |  `TransmissionProbabilityOnContact`               | Transmission risk for people located in the susceptible compartments. |
| $\xi_{I_{NS}}$               |  `RelativeTransmissionNoSymptoms`               | Proportion of nonsymptomatically infected people who are not isolated. |
| $\xi_{I_{Sy}}$               | `RiskOfInfectionFromSymptomatic`                | Proportion of infected persons with symptoms who are not isolated. |
| $N$                         | `m_N`   | Total population. |
| $D$                         |  Entry of `m_populations`  | Number of dead people. |
| $\mu_{z_1}^{z_2}$              |   `TransitionProbabilities`              | Probability of transition from a compartment $z_1$ to a compartment $z_2$. |  
| $\gamma_{z_1}^{z_2}(\tau)$              |   `TransitionDistributions`              | Expected proportion of people who are still in compartment $z_1$ $\tau$ days after entering this compartment and who will move to compartment $z_2$ later in the course of the disease. |  

The simulation runs in discrete time steps using a non-standard numerical scheme. This approach is based on the paper ["A non-standard numerical scheme for an age-of infection epidemic model" by Messina et al., Journal of Computational Dynamics, 2022](https://doi.org/10.3934/jcd.2021029). 

## Examples

An example can be found at:

- [IDE minimal example](../../examples/ide_secir.cpp)

## Initialization 

- The file [parameters_io](parameters_io.h) provides functionality to compute initial data for the IDE-SECIR model based on real data. An example for this initialization method can be found at [IDE initialization example](../../examples/ide_initialization.cpp).

- There are various options for initializing a fictional scenario. Regardless of the approach, you must provide a history of values for the transitions and additional information to compute the initial distribution of the population in the compartments. This information must be of the following type:  

    - You can state the number of total confirmed cases `total_confirmed_cases` at time $t_0$. The number of recovered people is set accordingly and the remaining values are derived.
    - You can set the number of people in the `Susceptible` compartment at time $t_0$ via `m_populations`. Initial values of the other compartments are derived.
    - You can set the number of people in the `Recovered` compartment at time $t_0$ via `m_populations`. Initial values of the other compartments are derived.
    - If none of the above is used, the force of infection formula and the values for the initial transitions are used as in [Messina et al (2022)](https://doi.org/10.3934/jcd.2021029) to set the `Susceptible`s. 
