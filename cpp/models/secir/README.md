# ODE-based / SEIR- and SECIR-type models

This module models and simulates the epidemic using an ODE-based SECIR-type model approach. Unlike the agent based model that uses an particular agents, this model simulates the spread of COVID19 in a population with subpopulations being in different compartments such as 'Exposed', 'Infected' or 'Recovered'.

## Structure

The model consists of the following major classes:
1. Populations: Generic class to create groups and/or compartments, i.e., subpopulations of a total population, using multiple dimensions. Dimensions can be age, infection state, income, gender, etc. 
2. Uncertain Matrix: Implements uncertain contact patterns between different (age) groups, one representation being a ContactMatrix(Group). This matrix contains contact patterns between different groups and `dampings` that model a change in contact patterns by a multiplicative factor at a given day.
3. Dampings: A `damping` object being the combination of a particular day and a multiplicative factor that changes the contact patterns. Dampings can be overwritten by or combined with dampings at later times. In order to not create discontinuities in pattern changes, the transition is smoothed by a cosine S-type function on an interval of maximum length of one day and which is C^1-smooth between two consecutive dampings.
5. SECIR: implements an *age-resolved ODE-model*, based on the non-age-resolved based model as described in https://www.medrxiv.org/content/10.1101/2020.04.04.20053637v2, uses the compartments 'Susceptible (S)', 'Exposed (E)', 'Carrier (C)', 'Infected (I)', 'Hospitalized (H)', 'ICU care (U)', 'Recovered (R)' and 'Dead'. Recovered people remain immune. Uses `populations` to model different 'groups' for a particular age-range (first dimension) and an infection state (second dimension). Parameters are set as SecirParams; they contain contact patterns in form of a UncertainContactMatrix and an extended set of epidemiologic parameters.
6. Parameter Space: Factory class for the SecirParams to set distributions to the different parameters and providing the opportunity to sample from these parameter set containing random distributions.
7. Parameter Studies: Method to be called on a set of SecirParams with a given set of random distributions to sample from the distributions and run ensemble run simulations with the obtained samples.

## Simulation

The simulation runs in discrete time steps using a numerical integration scheme. At each time step, a part of the population of each age-aware compartment moves from the current compartment to a new one. Different numerical integrations schemes are available, see the `math` folder. The SecirSimulation class handles the SecirParams and the numerical integrator. It also stores the result. Ensemble runs can be done using the Parameter Studies class as soon as random distributions are set for all the SecirParams. This can be done using the Parameter Space class.

## Examples

Different examples can be found at:

- examples/secir.cpp
- examples/secir_ageres.cpp
- examples/parameter_study_secir.cpp
