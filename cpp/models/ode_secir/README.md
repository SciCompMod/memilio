# ODE-based SECIR-type model

This module models and simulates the epidemic using an ODE-based SECIR-type model approach. Unlike the agent based model that uses particular agents, this model simulates the spread of a communicable disease in a population with subpopulations being in different compartments such as `Exposed`, `InfectedNoSymptoms`, `InfectedSymptoms` or `Recovered`. Individuals with and without symptoms are supposed to be infectious.

## Structure

The model consists of the following major classes:
1. Populations: Generic class to create groups and/or compartments, i.e., subpopulations of a total population, using multiple dimensions. Dimensions can be age, infection state, income, gender, etc. 
2. Uncertain Matrix: Implements uncertain contact patterns between different (age) groups, one representation being a ContactMatrix(Group). This matrix contains contact patterns between different groups and `Dampings` that model a change in contact patterns by a multiplicative factor at a given day.
3. Dampings: A `Damping` object is the combination of a particular day and a multiplicative factor that changes the contact patterns. Dampings can be overwritten by or combined with dampings at later times. In order to not create discontinuities in pattern changes, the transition is smoothed by a cosine S-type function on an interval of maximum length of one day. The resulting contact rate satisfies C^1-smoothness between two consecutive dampings.
4. SECIR: implements an *age-resolved ODE-model*, based on the non-age-resolved based model as described in https://www.medrxiv.org/content/10.1101/2020.04.04.20053637v1, uses the compartments `Susceptible (S)`, `Exposed (E)`, `InfectedNoSymptoms (I_NS)`, `InfectedSymptoms (I_Sy)`, `InfectedSevere (I_Sev)`, `InfectedCritical (I_Cr)`, `Recovered (R)` and `Dead`. The final model has been published in https://doi.org/10.1016/j.mbs.2021.108648. `Exposed` individuals are not infectious and it is assumed that subsequent infections can occur at random during the infectious period before the onset of symptoms, i.e., from state `InfectedNoSymptoms (I_NS)`. Recovered people remain immune. Severely or critically infected individuals are assumed to be isolated. The `Model` Uses `Populations` to model different 'groups' of a particular age-range (first dimension) and an `InfectionState` (second dimension). Parameters are set as `Parameters`; they contain contact patterns in form of an `UncertainContactMatrix` and an extended set of pathogen-dependent parameters.
5. Parameter Space: Factory class for the 'Parameters' to set distributions to the different parameters and providing the opportunity to sample from these parameter set containing random distributions.
6. Parameter Studies: Method to be called on a set of 'Parameters' with a given set of random distributions to sample from the distributions and run ensemble simulations with the obtained samples.

Below is an overview of the model architecture and its compartments.

![secir_model](https://github.com/SciCompMod/memilio/assets/70579874/46b09e8a-d083-4ef9-8328-21975890b60f)
| Mathematical variable                   | C++ variable name | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\phi$                      |  `ContactPatterns`               | Matrix of daily contact rates / number of daily contacts between different age groups. |
| $\rho$                      |  `TransmissionProbabilityOnContact`               | Transmission risk for people located in one of the susceptible compartments. |
| $\xi_{I_{NS}}$               |  `RelativeTransmissionNoSymptoms`               | Proportion of nonsymptomatically infected people who are not isolated. |
| $\xi_{I_{Sy}}$               | `riskFromInfectedSymptomatic`                | Proportion of infected people with symptomps who are not isolated (time-dependent if `TestAndTraceCapacity` used). |
| $N_j$                         | `Nj`   | Total population of age group $j$. |
| $D_i$                         |  `Di`  | Number of death people of age group $i$. |
| $\frac{1}{T_{E}}$                    |  `rateE`               | Time in days an individual stays in the Exposed compartment  (Computed from `SerialInterval` and `IncubationPeriod`). |
| $\frac{1}{T_{I_{NS}}}$                    |  `rateINS`               | Time in days an individual stays in the InfectedNoSymptoms compartment (Computed from `SerialInterval` and `IncubationPeriod`). |
| $T_{I_{Sy}}$                    |  `TimeInfectedSymptoms`               | Time in days an individual stays in the InfectedSymptoms compartment. |
| $T_{I_{Sev}}$                       |  `TimeInfectedSevere`               | Time in days an individual stays in the InfectedSevere compartment. |
| $T_{I_{Cr}}$                       |  `TimeInfectedCritical`               | Time in days an individual stays in the InfectedCritical compartment. |
| $\mu_{I_{NS}}^{I_{Sy}}$              |   `1 - RecoveredPerInfectedNoSymptoms`              | Probability of transition from compartment InfectedNoSymptoms to InfectedSymptoms. |  
| $\mu_{I_{Sy}}^{I_{Sev}}$              |   `SeverePerInfectedSymptoms`              | Probability of transition from compartment InfectedSymptoms to InfectedSevere. |
| $\mu_{I_{Sev}}^{I_{Cr}}$              |   `CriticalPerSevere`              | Probability of transition from compartment InfectedSevere to InfectedCritical. |  
| $\mu_{I_{Cr}}^{D}$              |   `DeathsPerCritical`              | Probability of dying when in compartment InfectedCritical. |   


## Simulation

The simulation runs in discrete time steps using a numerical integration scheme. At each time step, a part of the population of each age-aware compartment moves from the current compartment to a new one. Different numerical integrations schemes are available, see the `math` folder. The Simulation class handles the parameters and the numerical integrator. It also stores the result. Ensemble runs can be done using the Parameter Studies class as soon as random distributions are set for all the parameters. This can be done using the Parameter Space class.

## Examples

Different examples can be found at:

- examples/ode_secir.cpp
- examples/ode_secir_ageres.cpp
- examples/ode_secir_parameter_study.cpp
