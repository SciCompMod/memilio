# Stochastic Metapopulation Model

The stochastic metapopulation model uses a Markov process to simulate disease dynamics. Agents have an infection state and a position which is given by the region index the agent is in, hence agents in the same region have the same position. The evolution of the system state is determined by the following master equation

```math
\partial_t p(X,Z;t) = G p(X,Z;t) + L p(X,Z;t).
```
The operator $G$ defines the infection state adoptions and only acts on $Z$, while $L$ defines location changes, only acting on $X$. Infection state adoptions are modeled as stochastic jumps with independent Poisson processes given by adoption rate functions. As agents' positions are given by their subregion index, there is no movement within one subregion. Movement between subregions is also modeled as stochastic jumps with independent Poisson processes given by transition rate functions. Gillespie's direct methode (stochastic simulation algorithm) is used for simulation.

## Simulation

At the beginning of the simulation, the waiting times for all events are drawn. Then the time is advanced until the time point of the next event - which can be a spatial transition or an infection state adoption - and the event takes places. The waiting times of the other events are updated and a new waiting time for the event that just happend is drawn. The simulation saves the system state in discrete time steps.

For a detailed description and application of the model, see:

- Bicker J, Schmieding R, et al. (2025) Hybrid metapopulation agent-based epidemiological models for efficient insight on the individual scale: A contribution to green computing.
Infectious Disease Modelling, Volume 10, Issue 2. https://doi.org/10.1016/j.idm.2024.12.015
