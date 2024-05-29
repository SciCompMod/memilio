# GLCT SECIR model

This model is based on the Generalized Linear Chain Trick (GLCT). 

The GLCT provides the option to use Phase-type distributed stay times in the compartments through the use of subcompartments. The Generalized Linear Chain Trick is an extension of the Linear Chain Trick (as the name already suggests). Phase-type distributions are dense in the field of all positive-valued distributions. Therefore, for any positive-valued distribution, a phase-type distribution of arbitrary precision can be identified.
The normal ODE models have (possibly unrealistic) exponentially distributed stay times.
The GLCT model can still be described by an ordinary differential equation system.

For the concept see 
- P. J. Hurtado and A. S. Kirosingh, "Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models“, 2019. (https://doi.org/10.1007/s00285-019-01412-w)
- P. J. Hurtado and C. Richards, "Building mean field ODE models using the generalized linear chain trick & Markov chain theory", 2021. (https://doi.org/10.1080/17513758.2021.1912418)

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
It is possible to include phase-type distributed stay times for the five compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.
