# ODE-based SECIR-type model

This module models and simulates the epidemic using an ODE-based SECIR-type model approach. Unlike the agent based model that uses an particular agents, this model simulates the spread of COVID19 in a population with subpopulations being in different compartments such as 'Exposed', 'Infected' or 'Recovered'.

## Structure

The model consists of the following major classes:
1. Populations: Generic class to create groups and/or compartments, i.e., subpopulations of a total population, using multiple dimensions. Dimensions can be age, infection state, income, gender, etc. 
2. Uncertain Matrix: Implements uncertain contact patterns between different (age) groups, one representation being a ContactMatrix(Group). This matrix contains contact patterns between different groups and `dampings` that model a change in contact patterns by a multiplicative factor at a given day.
3. Dampings: A `damping` object being the combination of a particular day and a multiplicative factor that changes the contact patterns. Dampings can be overwritten by or combined with dampings at later times. In order to not create discontinuities in pattern changes, the transition is smoothed by a cosine S-type function on an interval of maximum length of one day and which is C^1-smooth between two consecutive dampings.
5. SECIR: implements an *age-resolved ODE-model*, based on the non-age-resolved based model as described in https://www.medrxiv.org/content/10.1101/2020.04.04.20053637v2, uses the compartments 'Susceptible (S)', 'Exposed (E)', 'InfectedNoSymptoms (I_NS)', 'InfectedSymptoms (I_Sy)', 'InfectedSevere (I_Sev)', 'InfectedCritical (I_Cr)', 'Recovered (R)' and 'Dead'. Recovered people remain immune. Uses `populations` to model different 'groups' for a particular age-range (first dimension) and an infection state (second dimension). Parameters are set as 'Parameters'; they contain contact patterns in form of a UncertainContactMatrix and an extended set of epidemiologic parameters.
6. Parameter Space: Factory class for the 'Parameters' to set distributions to the different parameters and providing the opportunity to sample from these parameter set containing random distributions.
7. Parameter Studies: Method to be called on a set of 'Parameters' with a given set of random distributions to sample from the distributions and run ensemble run simulations with the obtained samples.

Below is an overview of the model architecture and its compartments.

![secir_model](https://github.com/DLR-SC/memilio/assets/69154294/9e0ba2fb-f966-442f-86ca-b568d74b9311)
| Parameter                   | Implementation | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\phi$                      |  `ContactPatterns`               | Daily contact rate |
| $\rho$                      |  `TransmissionProbabilityOnContact`               | Transmission risk for people located in one of the susceptible compartments |
| $\xi_{I_{NS}}$               |  `RelativeTransmissionNoSymptoms`               | Proportion of asymptomatic infected people who are not isolated |
| $\xi_{I_{Sy}}$               | `riskFromInfectedSymptomatic`                | Proportion of infected people with symptomps who are not isolated. |
| $N$                         | `Nj`   | Total population. |
| $D_*$                         |  `-`  | Number of death people. `Nj` is directly calculated excluding death people. |
| $\frac{1}{T_{E}}$                    |  `rateE`               | Time in days an individual stays in the Exposed compartment. |
| $\frac{1}{T_{I_{NS}}}$                    |  `rateINS`               | Time in days an individual stays in the Infected No Symptoms compartment. |
| $T_{I_{Sy}}$                    |  `TimeInfectedSymptoms`               | Time in days an individual stays in the Infected Symptoms compartment. |
| $T_{I_{Sev}}$                       |  `TimeInfectedSevere`               | Time in days an individual stays in the Infected Severe compartment. |
| $T_{I_{Cr}}$                       |  `TimeInfectedCritical`               | Time in days an individual stays in the Infected Critical compartment. |
| $\mu_{I_{NS}}^{I_{Sy}}$              |   `1 - RecoveredPerInfectedNoSymptoms`              | Probability of transition from compartment Infected No Symptoms to Infected Symptoms. |  
| $\mu_{I_{Sy}}^{I_{Sev}}$              |   `SeverePerInfectedSymptoms`              | Probability of transition from compartment Infected Symptoms to Infected Severe. |
| $\mu_{I_{Sev}}^{I_{Cr}}$              |   `CriticalPerSevere`              | Probability of transition from compartment Infected Severe to Infected Critical. |  
| $\mu_{I_{Cr}}^{D}$              |   `DeathsPerCritical`              | Probability of dying when located in compartment Infected Critical. |   


## Simulation

The simulation runs in discrete time steps using a numerical integration scheme. At each time step, a part of the population of each age-aware compartment moves from the current compartment to a new one. Different numerical integrations schemes are available, see the `math` folder. The Simulation class handles the parameters and the numerical integrator. It also stores the result. Ensemble runs can be done using the Parameter Studies class as soon as random distributions are set for all the parameters. This can be done using the Parameter Space class.

## Examples

Different examples can be found at:

- examples/ode_secir.cpp
- examples/ode_secir_ageres.cpp
- examples/ode_secir_parameter_study.cpp
