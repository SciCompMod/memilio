# LCT SECIR TWO DISEASES model

This model describes infection with two independent (e.g. no co-infection) diseases a and b, based on the Linear Chain Trick. 

The Linear Chain Trick (LCT) provides the option to use Erlang-distributed stay times in the compartments through the use of subcompartments. 
The normal ODE models have (possibly unrealistic) exponentially distributed stay times.
The LCT model can still be described by an ordinary differential equation system.

For the concept of LCT models with one disease see 
- Lena Plötzke, "Der Linear Chain Trick in der epidemiologischen Modellierung als Kompromiss zwischen gewöhnlichen und Integro-Differentialgleichungen", 2023. (https://elib.dlr.de/200381/, German only)
- P. J. Hurtado und A. S. Kirosingh, "Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models“, 2019. (https://doi.org/10.1007/s00285-019-01412-w)

For each infection there are the infection states Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere, InfectedCritical, Recovered and Dead
The infectious compartments are InfectedNoSymptoms and InfectedSymptoms, all other compartments are considered to be not infectious

There are two possibilities for a susceptible individual (since we assume no co-infection): 
1. Get infected with disease a, then (if not dead) get infected with disease b
2. Get infected with disease b, then (if not dead) get infected with disease a

This leads to the following compartments:
- `Susceptible` ($S$), may get infected with a or b (-> Exposed_1a or Exposed_1b) at any time
Infection 1a:
- `Exposed_1a` ($E_{1a}$), becomes InfectedNoSymptoms_1a after some time
- `InfectedNoSymptoms_1a` ($I_{NS, 1a}_$), becomes InfectedSymptoms_1a or Recovered_a after some time
- `InfectedSymptoms_1a` ($I_{Sy, 1a}$), becomes InfectedSevere_1a or Recovered_a after some time
- `InfectedSevere_1a` ($I_{Sev, 1a}$), becomes InfectedCritical_1a or Recovered_a after some time
- `InfectedCritical_1a` ($I_{Cr, 1a}$), becomes Recovered_a or Dead_a after some time
- `Recovered_a` ($R_a$), immune to a, may get infected with b (-> Exposed_2b) at any time
- `Dead_a` ($D_a$), absorbing state
Infection 1b:
- `Exposed_1b` ($E_{1b}$), becomes InfectedNoSymptoms_1b after some time
- `InfectedNoSymptoms_1b` ($I_{NS, 1b}$), becomes InfectedSymptoms_1b or Recovered_1b after some time
- `InfectedSymptoms_1b` ($I_{Sy, 1b}$), becomes InfectedSevere_1b or Recovered_b after some time
- `InfectedSevere_1b` ($I_{Sev, 1b}$), becomes InfectedCritical_1b or Recovered_b after some time
- `InfectedCritical_1b` ($I_{Cr, 1b}$), becomes Recovered_b or Dead_b after some time
- `Recovered_b` ($R_b$), immune to b, may get infected with b (-> Exposed_2a) at any time
- `Dead_b` ($D_b$), absorbing state
Infection 2a:
- `Exposed_2a` ($E_{2a}$), becomes InfectedNoSympotoms_2a after some time
- `InfectedNoSymptoms_2a` ($I_{NS, 2a}$), becomes InfectedSymptoms_2a or Recovered_ab after some time
- `InfectedSymptoms_2a` ($I_{Sy, 2a}$), becomes InfectedSevere_2a or Recovered_ab after some time
- `InfectedSevere_2a` ($I_{Sev, 2a}$), becomes InfectedCritical_2a or Recovered_ab after some time
- `InfectedCritical_2a` ($I_{Cr, 2a}$), becomes Recovered_ab or Dead_a after some time
- `Recovered_ab` ($R_{ab}$), absorbing state, immune to a and b
Infection 2b:
- `Exposed_2b` ($E_{2b}$), becomes infected after some time
- `InfectedNoSymptoms_2b` ($I_{NS,2b}$), becomes InfectedSymptoms_2b or Recovered_ab after some time
- `InfectedSymptoms_2b` ($I_{Sy,2b}$), becomes InfectedSevere_2b or Recovered_ab after some time
- `InfectedSevere_2b` ($I_{Sev, 2b}$), becomes InfectedCritical_2b or Recovered_ab after some time
- `InfectedCritical_2b` ($I_{Cr, 2b}$), becomes Recovered_ab or Dead_b after some time


It is possible to include subcompartments for the five compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.
You can divide the population according to different groups, e.g. AgeGroups or gender and choose parameters according to groups.

The parameters depend on the disease (a or b), the number of subcompartments depends on the individual compartment (since it can be set independently).


| Mathematical variable       | C++ variable name                      | Description                                                                                 |
|---------------------------- | -------------------------------------- | ------------------------------------------------------------------------------------------- |
| $\phi$                      |  `ContactPatterns`                     | Average number of contacts of a person per day.                                             |
| $\rho_i$                    |  `TransmissionProbabilityOnContact_i`  | Transmission risk for people located in the susceptible compartments, i in {a,b}.           |
| $\xi_{I_{NS},i}$            |  `RelativeTransmissionNoSymptoms_i`    | Proportion of nonsymptomatically infected people who are not isolated, i in {a,b}.          |
| $\xi_{I_{Sy},i}$            |  `RiskOfInfectionFromSymptomatic_i`    | Proportion of infected people with symptoms who are not isolated, i in {a,b}.               |
| $n_{E,i}$                   |   Defined in `LctStates`               | Number of subcompartments of the Exposed compartment, i in {1a, 2a, 1b,2b}.                 |
| $n_{NS,i}$                  |   Defined in `LctStates`               | Number of subcompartments of the InfectedNoSymptoms compartment, i in {1a, 2a, 1b,2b}.      |
| $n_{Sy,i}$                  |   Defined in `LctStates`               | Number of subcompartments of the InfectedSymptoms compartment, i in {1a, 2a, 1b,2b}.        |
| $n_{Sev,i}$                 |   Defined in `LctStates`               | Number of subcompartments of the InfectedSevere compartment, i in {1a, 2a, 1b,2b}.          |
| $n_{Cr,i}$                  |   Defined in `LctStates`               | Number of subcompartments of the InfectedCritical compartment, i in {1a, 2a, 1b,2b}.        |
| $T_{E,i}$                   |  `TimeExposed`                         | Average time in days an individual stays in the Exposed compartment, i in {a,b}.            |
| $T_{I_{NS},i}$              |  `TimeInfectedNoSymptoms`              | Average time in days an individual stays in the InfectedNoSymptoms compartment, i in {a,b}. |
| $T_{I_{Sy},i}$              |  `TimeInfectedSymptoms`                | Average time  in days an individual stays in the InfectedSymptoms compartmen, i in {a,b}t.  |
| $T_{I_{Sev},i}$             |  `TimeInfectedSevere`                  | Average time in days an individual stays in the InfectedSevere compartment, i in {a,b}.     |
| $T_{I_{Cr},i}$              |  `TimeInfectedCritical`                | Average time in days an individual stays in the InfectedCritical compartment, i in {a,b}.   |
| $\mu_{I_{NS},i}^{R}$        |  `RecoveredPerInfectedNoSymptoms`      | Probability of transition from compartment InfectedNoSymptoms to Recovered, i in {a,b}.     |  
| $\mu_{I_{Sy},i}^{I_{Sev}}$  |  `SeverePerInfectedSymptoms`           | Probability of transition from compartment InfectedSymptoms to InfectedSevere, i in {a,b}.  |
| $\mu_{I_{Sev},i}^{I_{Cr}}$  |  `CriticalPerSevere`                   | Probability of transition from compartment InfectedSevere to InfectedCritical, i in {a,b}.  |  
| $\mu_{I_{Cr},i}^{D}$        |  `DeathsPerCritical`                   | Probability of dying when in compartment InfectedCritical, i in {a,b}.                      |   


## Examples

A simple example can be found at [LCT2D minimal example](../../examples/lct_secir_2_diseases.cpp).
