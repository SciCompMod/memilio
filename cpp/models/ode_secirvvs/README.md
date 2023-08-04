# SECIR model with COVID-19 variants and vaccinations

This model extends the basic SECIR model by adding vaccinations and allowing the implicit modeling of a newly arriving variant that takes hold.

Vaccinations are modeled by adding compartments for partially and fully vaccinated persons. `Partially and fully vaccinated` is to be understood in this context as the person having received a first and second vaccine shot as in 2021. These model lines can be reused by resetting parameters. Persons that have recovered from the disease are treated as fully vaccinated from that time forward. Vaccinated persons are added on every day of simulation, see parameters `DailyFirstVaccination` and `DailyFullVaccination`. All groups can get an infection or get reinfected. Vaccinated persons are less likely to develop symptoms. E.g., the probability to develop symptoms when carrying the virus is the base probability from the SECIR model multiplied with the `ReducInfectedSymptomsPartialImmunity` parameter.

The ratio of two variants can change over time, which affects the average transmissiblity of the disease. Infectiousness of different variants can be set in the parameters.

## Examples

The extended model is used in the 2021_vaccination_sarscov2_delta_germany simulation. 
An easier example can be found in [examples/ode_secirvvs.cpp](../../examples/ode_secirvvs.cpp)

Examples of the basic SECIR model can be found at:

- examples/ode_secir.cpp
- examples/ode_secir_ageres.cpp
- examples/ode_secir_parameter_study.cpp
