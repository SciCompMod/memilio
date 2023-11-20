# IDE SECIR model

This model is based on Integro-differential equations. The eight compartments 
- Susceptible, may become exposed at any time
- Exposed, becomes infected after some time
- InfectedNoSymptoms, becomes InfectedSymptoms or Recovered after some time
- InfectedSymptoms, becomes InfectedSevere or Recovered after some time
- InfectedSevere, becomes InfectedCritical or Recovered after some time
- InfectedCritical, becomes Recovered or Dead after some time
- Recovered
- Dead

are used to simulate the spread of the disease.


The simulation runs in discrete time steps using a non-standard numerical scheme. This approach is based on the paper "A non-standard numerical scheme for an age-of infection epidemic model" by Messina et al., Journal of Computational Dynamics, 2022. 

## Examples

An example can be found at:

- examples/ide_secir.cpp