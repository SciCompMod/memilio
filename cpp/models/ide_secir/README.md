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

| Mathematical variable                   | C++ variable name | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\phi$                      |  `ContactPatterns`               | Average number of contacts of a person per day. |
| $\rho$                      |  `TransmissionProbabilityOnContact`               | Transmission risk for people located in the susceptible compartments. |
| $\xi_{I_{NS}}$               |  `RelativeTransmissionNoSymptoms`               | Proportion of nonsymptomatically infected people who are not isolated. |
| $\xi_{I_{Sy}}$               | `RiskOfInfectionFromSymptomatic`                | Proportion of infected people with symptomps who are not isolated. |
| $N$                         | `m_N`   | Total population. |
| $D$                         |  `D`  | Number of dead people. |
| $\mu_{z_1}^{z_2}$              |   `TransitionProbabilities `              | Probability of transition from a compartment $z_1$ to a compartment $z_2$. |  
| $\gamma_{z_1}^{z_2}(\tau)$              |   `TransitionDistributions `              | TODO: The share of people that already spent $\tau$ days in compartment $z_1$ if they go to  $z_2$ in their course of infection. Should be a ember of StateAgeFunction.|  



The simulation runs in discrete time steps using a non-standard numerical scheme. This approach is based on the paper "A non-standard numerical scheme for an age-of infection epidemic model" by Messina et al., Journal of Computational Dynamics, 2022. 

## Examples

An example can be found at:

- [IDE minimal example](../../examples/ide_secir.cpp)

## Initialization 

- The file [parameters_io](parameters_io.h) provides functionality to compute initial data for the IDE-SECIR model based on real data. An example for this initialization method can be found at [IDE initialization example](../../examples/ide_initialization.cpp).