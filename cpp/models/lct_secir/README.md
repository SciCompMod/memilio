# LCT SECIR model

This model is based on the Linear Chain Trick. 
For the concept see 
- Lena Plötzke, "Der Linear Chain Trick in der epidemiologischen Modellierung als Kompromiss zwischen gewöhnlichen und Integro-Differentialgleichungen", 2023. (https://elib.dlr.de/200381/, German only)
- P. J. Hurtado und A. S. Kirosingh, "Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models“, 2019. (https://doi.org/10.1007/s00285-019-01412-w)

The eight compartments 
- Susceptible, may become exposed at any time
- Exposed, becomes infected after some time
- InfectedNoSymptoms, becomes InfectedSymptoms or Recovered after some time
- InfectedSymptoms, becomes InfectedSevere or Recovered after some time
- InfectedSevere, becomes InfectedCritical or Recovered after some time
- InfectedCritical, becomes Recovered or Dead after some time
- Recovered
- Dead

are used to simulate the spread of the disease. 
It is possible to include subcompartments for the five compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.

![tikzLCTSECIR](https://github.com/SciCompMod/memilio/assets/70579874/133bcec2-499c-4e95-af96-39299bf0fd67)

A simple example can be found at: examples/lct_secir.cpp. 


