LCT Models
==========

In MEmilio, two models based on Linear Chain Trick are implemented. Their respective usage is described below.


LCT-SECIR model
----------------

The Linear Chain Trick provides the option to use Erlang-distributed stay times in the compartments through the use of subcompartments. 
The normal ODE models have (possibly unrealistic) exponentially distributed stay times.
The LCT model can still be described by an ordinary differential equation system.

For the concept see 
- Lena Plötzke, "Der Linear Chain Trick in der epidemiologischen Modellierung als Kompromiss zwischen gewöhnlichen und Integro-Differentialgleichungen", 2023. (https://elib.dlr.de/200381/, German only)
- P. J. Hurtado und A. S. Kirosingh, "Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models“, 2019. (https://doi.org/10.1007/s00285-019-01412-w)

The eight compartments 
- `Susceptible` (:math:`S`), may become Exposed at any time
- `Exposed` (:math:`E`), becomes InfectedNoSymptoms after some time
- `InfectedNoSymptoms` (:math:`I_{NS}`), becomes InfectedSymptoms or Recovered after some time
- `InfectedSymptoms` (:math:`I_{Sy}`), becomes InfectedSevere or Recovered after some time
- `InfectedSevere` (:math:`I_{Sev}`), becomes InfectedCritical or Recovered after some time
- `InfectedCritical` (:math:`I_{Cr}`), becomes Recovered or Dead after some time
- `Recovered` (:math:`R`)
- `Dead` (:math:`D`)

are used to simulate the spread of the disease. 
It is possible to include subcompartments for the five compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.
You can divide the population according to different groups, e.g. AgeGroups or gender and choose parameters according to groups.


How to: Set up and run a simulation of the LCT-SECIR model
-----------------------------------------------------------




GLCT-SECIR model
-----------------


How to: Set up and run a simulation of the GLCT-SECIR model
------------------------------------------------------------
