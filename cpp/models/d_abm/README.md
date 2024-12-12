# Diffusive Agent-Based Model

This agent-based model uses a Markov process to simulate disease dynamics. The features of an agent are position and infection state. The evolution of the system state is determined by the following master equation

```math
\partial_t p(X,Z;t) = G p(X,Z;t) + L p(X,Z;t)
```
The operator $G$ defines the infection state adoptions and only acts on $Z$, while $L$ defines location changes, only acting on $X$. Infection state adoptions are modeled with independent Poisson processes given by adoption rate functions. Movement is modeled with independent diffusion processes. A temporal Gillespie algorithm is used for simulation.

The Model class needs an Implementation class as template argument which prvides the domain agents move and interact in. We here implemented a quadwell potential given in the class QuadWellModel, but any other suitable potential can be used as implementation. 

## Simulation

The simulation runs in discrete time steps. In every step we advance the model until the next infection state adoption event, then adopt the corresponding agent's infection state and draw a new waiting time until the next adoption event. If the waiting time until the next adoption event is bigger than the remaining time in the time step, we advance the model until the end of the time step.
