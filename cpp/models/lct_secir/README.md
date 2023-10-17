# LCT SECIR model

This model is based on the Linear Chain Trick. The eight compartments 
- Susceptible, may become exposed at any time
- Exposed, becomes infected after some time
- InfectedNoSymptoms, becomes InfectedSymptoms or Recovered after some time
- InfectedSymptoms, becomes InfectedSevere or Recovered after some time
- InfectedSevere, becomes InfectedCritical or Recovered after some time
- InfectedCritical, becomes Recovered or Dead after some time
- Recovered
- Dead

are used to simulate the spread of the disease. 
It ist possible to include subcompartments for the five compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.


## Examples

A simple example can be found at: examples/lct_secir.cpp. 
More complex examples with simulations of different scenarios can be found at: 
- examples/lct_secir_initializations.cpp
- examples/lct_secir_fictional_scenario.cpp 
- examples/lct_secir_fictional_scenario.cpp
The parameters used for simulations are calculated via examples/compute_parameters.cpp.
