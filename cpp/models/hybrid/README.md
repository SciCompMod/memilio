# Hybrid models

This model is a hybrid model combining two models of different types. The motivation behind combining different model types is that multiple modeling approaches have different advantages and disadvantages which can lead to a trade-off when chosing an appropriate model. For example, agent-based models offer a high level of detail while population-based models (e.g. based on ordinary differential equations) only deliver aggregated results. However, agent-based models need lots of data to be initialized and their runtime scales with the number of agents. Population-based models are in contrast computationally efficient even for large populations. To combine the advantages of different models, several hybridization approaches can be used.

## Temporal Hybridization

The temporal-hybrid model switches between two models during the course of the simulation according to a given condition. The switching condition can in general be chosen arbitrarily e.g. a fixed time point or a case number threshold can be used for the switch. What condition suits best depends on the used models and the concrete application. Both models have to be initialized with their corresponding parameters and are handed to the class TemporalHybridSimulation.

The TemporalHybridSimulation class needs the used model types as well as their result types as template arguments. The results of the models are used to evaluate the switching condition. Additionally, conversion functions to convert the first model to the second model and vice versa have to be implemented for the used model types.
We implemented conversion function for the following model combinations:

- Diffusive agent-based model, see mio::dabm, using the singlewell potential and the stochastic metapopulation model, see mio::smm.
- Diffusive agent-based model, see mio::dabm, using the singlewell potential, and the ODE-based SECIR-type model, see mio::osecir.

### Simulation

The simulation runs in discrete time steps. In every step, we first check whether we need to switch the model by evaluating the switching condition. If the condition is fulfilled, we convert the currently used model into the target model i.e. we use its current state/result to update the state of the target model. Then we advance the (new) currently used model until the next time step.

For a detailed description and application of the model, see:

- Bicker J, Schmieding R, et al. (2025) Hybrid metapopulation agent-based epidemiological models for efficient insight on the individual scale: A contribution to green computing.
Infectious Disease Modelling, Volume 10, Issue 2. https://doi.org/10.1016/j.idm.2024.12.015

