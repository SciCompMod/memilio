# LCT SECIR model

This model is based on the Linear Chain Trick. 

The Linear Chain Trick provides the Option to use Erlang-distributed sojourn times in the compartments through the use of subcompartments. 
The normal ODE models have (possibly unrealistic) exponentially distributed sojourn times.
The LCT model can still be described by an ordinary differential equation system.

For the concept see 
- Lena Plötzke, "Der Linear Chain Trick in der epidemiologischen Modellierung als Kompromiss zwischen gewöhnlichen und Integro-Differentialgleichungen", 2023. (https://elib.dlr.de/200381/, German only)
- P. J. Hurtado und A. S. Kirosingh, "Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models“, 2019. (https://doi.org/10.1007/s00285-019-01412-w)

The eight compartments 
- `Susceptible` (S), may become exposed at any time
- `Exposed` (E), becomes infected after some time
- `InfectedNoSymptoms` (C), becomes InfectedSymptoms or Recovered after some time
- `InfectedSymptoms` (I), becomes InfectedSevere or Recovered after some time
- `InfectedSevere` (H), becomes InfectedCritical or Recovered after some time
- `InfectedCritical` (U), becomes Recovered or Dead after some time
- `Recovered` (R)
- `Dead` (D)

are used to simulate the spread of the disease. 
It is possible to include subcompartments for the five compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.

Below is an overview of the model architecture and its compartments.

![tikzLCTSECIR](https://github.com/SciCompMod/memilio/assets/70579874/e385f26c-5070-4655-9faf-eba753dc8e99)

| Mathematical variable                   | C++ variable name | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\phi$                      |  `ContactPatterns`               | Average number of contacts of a person per day. |
| $\rho$                      |  `TransmissionProbabilityOnContact`               | Transmission risk for people located in the susceptible compartments. |
| $\xi_{C}$               |  `RelativeTransmissionNoSymptoms`               | Proportion of nonsymptomatically infected people who are not isolated. |
| $\xi_{I}$               | `RiskOfInfectionFromSymptomatic`                | Proportion of infected people with symptomps who are not isolated. |
| $N$                         | `m_N0`   | Total population. |
| $D$                         |  `D`  | Number of death people. |
| $n_E$                         |  Defined via `InfectionState`  | Number of subcompartments of the Exposed compartment. |
| $n_C$                         |  Defined via `InfectionState`  | Number of subcompartments of the InfectedNoSymptoms compartment. |
| $n_I$                         |  Defined via `InfectionState`  | Number of subcompartments of the InfectedSymptoms compartment. |
| $n_H$                         | Defined via `InfectionState`  | Number of subcompartments of the InfectedSevere compartment.|
| $n_U$                         |  Defined via `InfectionState`  | Number of subcompartments of the InfectedCritical compartment. |
| $T_E$                    |  `TimeExposed`               | Average time in days an individual stays in the Exposed compartment. |
| $T_C$                    |  `TimeInfectedNoSymptoms`               | Average time in days an individual stays in the InfectedNoSymptoms compartment. |
| $T_{I}$                    |  `TimeInfectedSymptoms`               | Average time  in days an individual stays in the InfectedSymptoms compartment. |
| $T_{H}$                       |  `TimeInfectedSevere`               | Average time in days an individual stays in the InfectedSevere compartment. |
| $T_{U}$                       |  `TimeInfectedCritical`               | Average time in days an individual stays in the InfectedCritical compartment. |
| $\mu_{C}^{R}$              |   `RecoveredPerInfectedNoSymptoms`              | Probability of transition from compartment InfectedNoSymptoms to Recovered. |  
| $\mu_{I}^{H}$              |   `SeverePerInfectedSymptoms`              | Probability of transition from compartment InfectedSymptoms to InfectedSevere. |
| $\mu_{H}^{U}$              |   `CriticalPerSevere`              | Probability of transition from compartment InfectedSevere to InfectedCritical. |  
| $\mu_{U}^{D}$              |   `DeathsPerCritical`              | Probability of dying when in compartment InfectedCritical. |   

The notation of the compartments with indices here stands for subcompartments and not for age groups or similar, as in some ODE models. Accordingly, $E_{n_E}$, for example, stands for the number of people in the $n_E$-th subcompartment of the Exposed compartment.


## Examples

A simple example can be found at [LCT minimal example](../../examples/lct_secir.cpp).

