# SECIRS-type model including multi-layer waning immunity

This model extends the SECIRVVS model by adding waning immunity and introducing temporary immunity states that change the meaning of recovery.

In the model, waning immunity is defined by the parameters `TimeWaningPartialImmunity`, `TimeWaningImprovedImmunity`, `TimeTemporaryImmunityPI`, and `TimeTemporaryImmunityII`. `TimeWaningPartialImmunity` and `TimeWaningImprovedImmunity` represent the (mean) duration after which an individual transitions  from one immunity layer to the next less protected one due to waning immunity, assuming no vaccination or recovery from infection has occurred during this period. The parameters `TimeTemporaryImmunityPI` and `TimeTemporaryImmunityII` denote the (mean) duration of temporary immunity following exposure to the virus, either through vaccination or recovery. During this state of temporary immunity, individuals are protected from reinfection and are incapable of transmitting the virus to others. Should individuals previously reside in the naive or partial immunity layer, their stay in the temporary immunity state results in a transition to the next more protected immunity layer.

For more details about the model, we refer to [1](https://www.medrxiv.org/content/10.1101/2024.03.01.24303602v1).

Below is an overview of the model architecture and its compartments.

![SECIRVVS_model](https://github.com/SciCompMod/memilio/assets/69154294/6dec331f-bd91-410f-be5e-c8cf6eb0572b)
| Mathematical variable                   | C++ variable name | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\lambda_{N,i} =  \rho_{N,i} \sum_j \phi_{i,j}\frac{\xi_{I_{NS}} (I_{NS,N,j} + I_{NS,PI,j} + I_{NS,II,j}) + \xi_{I_{Sy}} (I_{Sy,N,j} + I_{Sy,PI,j}+ I_{Sy,II,j})}{N_j^{D^\perp}}$                      |  `ext_inf_force_dummy`               | Force of infection for susceptibles located in the naive immunity level. |
| $\lambda_{PI,i} = \rho_{PI,i}\sum_j \phi_{i,j}\frac{\xi_{I_{NS}} (I_{NS,N,j} + I_{NS,PI,j} + I_{NS,II,j}) + \xi_{I_{Sy}} (I_{Sy,N,j} + I_{Sy,PI,j}+ I_{Sy,II,j})}{N_j^{D^\perp}}$                      |  `reducExposedPartialImmunity * ext_inf_force_dummy`               | Force of infection for susceptibles located in the partial immunity level. |
| $\lambda_{II,i} = \rho_{II}\sum_j \phi_{i,j}\frac{\xi_{I_{NS}} (I_{NS,N,j} + I_{NS,PI,j} + I_{NS,II,j}) + \xi_{I_{Sy}} (I_{Sy,N,j} + I_{Sy,PI,j}+ I_{Sy,II,j})}{N_j^{D^\perp}}$                     |  `reducExposedImprovedImmunity * ext_inf_force_dummy`               | Force of infection for susceptibles located in the improved immunity level. |
| $\phi$                      |  `ContactPatterns`               | Matrix of daily contact rates / number of daily contacts between different age groups. |
| $\rho$                      |  `TransmissionProbabilityOnContact`               | Transmission risk for people located in one of the susceptible compartments. |
| $\xi_{I_{NS}}$               |  `RelativeTransmissionNoSymptoms`               | Proportion of asymptomatically infected people who are not isolated (time-dependent if `TestAndTraceCapacity` used). |
| $\xi_{I_{Sy}}$               | `riskFromInfectedSymptomatic`                | Proportion of symptomatically infected people who are not isolated (time-dependent if `TestAndTraceCapacity` used). |
| $N_j^{D^\perp}$                         | `Nj`   | Sum of all living individuals of age groups j. |
| $T_{E}$                    |  `TimeExposed`               | Time in days an individual stays in the Exposed compartment. |
| $T_{I_{NS}}$                    |  `TimeInfectedNoSymptoms`               | Time in days an individual stays in the InfectedNoSymptoms compartment. |
| $T_{I_{Sy}}$                    |  `TimeInfectedSymptoms`               | Time in days an individual stays in the InfectedSymptoms compartment. |
| $T_{I_{Sev}}$                       |  `TimeInfectedSevere`               | Time in days an individual stays in the InfectedSevere compartment. |
| $T_{I_{Cr}}$                       |  `TimeInfectedCritical`               | Time in days an individual stays in the InfectedCritical compartment. |
| $T_{\mathcal{I}_{PI}}$                       |  `TimeTemporaryImmunityPI`               | Time in days an individual stays in the TemporaryImmunPartialImmunity compartment. |
| $T_{\mathcal{I}_{PI}}$                       |  `TimeTemporaryImmunityII`               | Time in days an individual stays in the TemporaryImmunImprovedImmunity compartment. |
| $T_{W_{PI}}$                       |  `TimeWaningPartialImmunity`               | Time in days an individual stays in the SusceptiblePartialImmunity compartment before waning to the SusceptibleNaive compartment assuming no exposure occurred during this period. |
| $T_{W_{II}}$                       |  `TimeWaningImprovedImmunity`               | Time in days an individual stays in the SusceptibleImprovedImmunity compartment before waning to the SusceptiblePartialImmunity compartment assuming no exposure occurred during this period. |
| $\mu_{I_{NS}}^{I_{Sy}}$              |   `1 - RecoveredPerInfectedNoSymptoms`              | Probability of transition from compartment InfectedNoSymptoms to InfectedSymptoms. |  
| $\mu_{I_{Sy}}^{I_{Sev}}$              |   `SeverePerInfectedSymptoms`              | Probability of transition from compartment InfectedSymptoms to InfectedSevere. |
| $\mu_{I_{Sev}}^{I_{Cr}}$              |   `CriticalPerSevere`              | Probability of transition from compartment InfectedSevere to InfectedCritical. |  
| $\mu_{I_{Cr}}^{D}$              |   `DeathsPerCritical`              | Probability of dying when located in compartment InfectedCritical. |   
| $\kappa$              |   `ReducTimeInfectedMild`              | Reduction factor for time intervals for specific partial and improved immunity compartments. |   