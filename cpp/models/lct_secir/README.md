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
It ist possible to include subcompartments for the five compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere,and InfectedCritical.


## Examples

An example can be found at:
