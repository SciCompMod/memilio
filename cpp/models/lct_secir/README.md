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

A simple example can be found at: examples/lct_secir.cpp. 

The file [initializer_flows](initializer_flows.h) provides functionality to compute an initial value vector for the LCT-SECIR model based on initial data in the form of a TimeSeries of InfectionTransition%s. For the concept of the InfectionTransition%s or flows, see also the IDE-SECIR model. This method can be particularly useful if a comparison is to be made with an IDE model with matching initialization or if the real data is in the form of flows.


