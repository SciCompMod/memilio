GLCT model
===========

Introduction
-------------

This model is based on the Generalized Linear Chain Trick (GLCT). 

The GLCT provides the option to use phase-type distributed stay times in the compartments through the use of subcompartments. The Generalized Linear Chain Trick is an extension of the Linear Chain Trick (as the name already suggests). Phase-type distributions are dense in the field of all positive-valued distributions. Therefore, for any positive-valued distribution, a phase-type distribution of arbitrary precision can be identified.
The normal ODE models have (possibly unrealistic) exponentially distributed stay times.
The GLCT model can still be described by an ordinary differential equation system.

For the concept see 

- Hurtado PJ and Kirosingh AS (2019) Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models. Journal of Mathematical Biology. https://doi.org/10.1007/s00285-019-01412-w
- Hurtado PF and Richards C (2021) Building mean field ODE models using the generalized linear chain trick & Markov chain theory. Journal of Biological Dynamics. https://doi.org/10.1080/17513758.2021.1912418

Here, the eight compartments 

- `Susceptible` (:math:`S`), may become Exposed at any time
- `Exposed` (:math:`E`), becomes InfectedNoSymptoms after some time
- `InfectedNoSymptoms` (:math:`I_{NS}`), becomes InfectedSymptoms or Recovered after some time
- `InfectedSymptoms` (:math:`I_{Sy}`), becomes InfectedSevere or Recovered after some time
- `InfectedSevere` (:math:`I_{Sev}`), becomes InfectedCritical or Recovered after some time
- `InfectedCritical` (:math:`I_{Cr}`), becomes Recovered or Dead after some time
- `Recovered` (:math:`R`)
- `Dead` (:math:`D`)

are used to simulate the spread of the disease. 
It is possible to include phase-type distributed stay times for the five compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.

Simulation
-----------

How to: Set up and run a simulation of the GLCT-SECIR model
------------------------------------------------------------
