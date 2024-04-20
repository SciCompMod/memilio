
# SDE SEIR two variant compartment model

This model is a stochastic multivariant SDE SEIR model with  10 compartments and few parameters, addressing reinfection with a second variant. 
- Susceptible, may become exposed at any time.
- Exposed, may become infectious at any time
- Infected, will be recovered after some time.
- Recovered, recovered from previous infection, 

We assume simulations over short periods of time, so that the population size can be considered constant and birth as well as (natural) mortality rates can be ignored. 

Important note on the solution: The numerical integration method is a Euler-Maruyama which reuses the EulerIntegratorCore of MEmilio accordingly. The (re)use of other implemented integrators for ODEs is neither considered nor suggested at the moment.



