# SECIR model with COVID-19 variants and vaccinations

This model extends the basic SECIR model by adding vaccinations and allowing the implicit modeling of a newly arriving variant that takes hold.

Vaccinations are modeled by adding compartments for partially and fully vaccinated persons. `Partially and fully vaccinated` is to be understood in this context as the person having received a first and second vaccine shot as in 2021. These model lines can be reused by resetting parameters. Persons that have recovered from the disease are treated as fully vaccinated from that time forward. Vaccinated persons are added on every day of simulation, see parameters `DailyFirstVaccination` and `DailyFullVaccination`. All groups can get an infection or get reinfected. Vaccinated persons are less likely to develop symptoms. E.g., the probability to develop symptoms when carrying the virus is the base probability from the SECIR model multiplied with the `ReducInfectedSymptomsPartialImmunity` parameter.

The ratio of two variants can change over time, which affects the average transmissiblity of the disease. Infectiousness of different variants can be set in the parameters.

Below is an overview of the model architecture and its compartments.

![SECIRVVS_model](https://github.com/DLR-SC/memilio/assets/69154294/cf5ffd74-245c-4558-9de2-71b82cf79441)
| Parameter                   | Implementation | Description |
|---------------------------- | --------------- | -------------------------------------------------------------------------------------------------- |
| $\lambda_{N,i} =  \rho_{N,i} \sum_j \phi_{i,j}\frac{\xi_{I_{NS}} (I_{NS,N,j} + I_{NS,PI,j}+ I_{NS,II,j}) + \xi_{I_{Sy}} (I_{Sy,N,j} + I_{Sy,PI,j}+ I_{Sy,II,j})}{N_j^{D^\perp}}$                      |  `ext_inf_force_dummy`               | Force of infection for susceptibles located in the naive immunity level |
| $\lambda_{PI,i} = \rho_{PI,i}\sum_j \phi_{i,j}\frac{\xi_{I_{NS}} (I_{NS,N,j} + I_{NS,PI,j}+ I_{NS,II,j}) + \xi_{I_{Sy}} (I_{Sy,N,j} + I_{Sy,PI,j}+ I_{Sy,II,j})}{N_j^{D^\perp}}$                      |  `reducExposedPartialImmunity * ext_inf_force_dummy`               | Force of infection for susceptibles located in the partial immunity level |
| $\lambda_{II,i} = \rho_{II}\sum_j \phi_{i,j}\frac{\xi_{I_{NS}} (I_{NS,N,j} + I_{NS,PI,j}+ I_{NS,II,j}) + \xi_{I_{Sy}} (I_{Sy,N,j} + I_{Sy,PI,j}+ I_{Sy,II,j})}{N_j^{D^\perp}}$                     |  `reducExposedImprovedImmunity * ext_inf_force_dummy`               | Force of infection for susceptibles located in the improved immunity level |
| $\phi_{i,j}$                      |  `ContactPatterns`               | Daily contact rate between two age groups i and j. |
| $\rho$                      |  `TransmissionProbabilityOnContact`               | Transmission risk for people located in one of the susceptible compartments |
| $\xi_{I_{NS}}$               |  `RelativeTransmissionNoSymptoms`               | Proportion of asymptomatic infected people who are not isolated |
| $\xi_{I_{Sy}}$               | `riskFromInfectedSymptomatic`                | Proportion of infected people with symptomps who are not isolated. |
| $N_j^{D^\perp}$                         | `Nj`   | Sum of all living individuals of age groups j. |
| $\frac{1}{T_{E}}$                    |  `rateE`               | Time in days an individual stays in the Exposed compartment. |
| $\frac{1}{T_{I_{NS}}}$                    |  `rateINS`               | Time in days an individual stays in the Infected No Symptoms compartment. |
| $T_{I_{Sy}}$                    |  `TimeInfectedSymptoms`               | Time in days an individual stays in the Infected Symptoms compartment. |
| $T_{I_{Sev}}$                       |  `TimeInfectedSevere`               | Time in days an individual stays in the Infected Severe compartment. |
| $T_{I_{Cr}}$                       |  `TimeInfectedCritical`               | Time in days an individual stays in the Infected Critical compartment. |
| $\mu_{I_{NS}}^{I_{Sy}}$              |   `1 - RecoveredPerInfectedNoSymptoms`              | Probability of transition from compartment Infected No Symptoms to Infected Symptoms. |  
| $\mu_{I_{Sy}}^{I_{Sev}}$              |   `SeverePerInfectedSymptoms`              | Probability of transition from compartment Infected Symptoms to Infected Severe. |
| $\mu_{I_{Sev}}^{I_{Cr}}$              |   `CriticalPerSevere`              | Probability of transition from compartment Infected Severe to Infected Critical. |  
| $\mu_{I_{Cr}}^{D}$              |   `DeathsPerCritical`              | Probability of dying when located in compartment Infected Critical. |   
| $\kappa$              |   `ReducTimeInfectedMild`              | Reduction factor for time intervals for specific partial and improved immunity compartments. |   

## Examples

The extended model is used in the 2021_vaccination_sarscov2_delta_germany simulation. 
An easier example can be found in [examples/ode_secirvvs.cpp](../../examples/ode_secirvvs.cpp)

Examples of the basic SECIR model can be found at:

- examples/ode_secir.cpp
- examples/ode_secir_ageres.cpp
- examples/ode_secir_parameter_study.cpp