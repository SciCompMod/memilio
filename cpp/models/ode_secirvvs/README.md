# SECIR model with COVID-19 variants and vaccinations

This model extends the basic SECIR model by adding vaccinations and allowing the implicit modeling of a newly arriving variant that takes hold.

Vaccinations are modeled by adding compartments for partially and fully vaccinated persons. `Partially and fully vaccinated` is to be understood in this context as the person having received a first and second vaccine shot as in 2021. These model lines can be reused by resetting parameters. Persons that have recovered from the disease are treated as fully vaccinated from that time forward. Vaccinated persons are added on every day of simulation, see parameters `DailyFirstVaccination` and `DailyFullVaccination`. All groups can get an infection or get reinfected. Vaccinated persons are less likely to develop symptoms. E.g., the probability to develop symptoms when carrying the virus is the base probability from the SECIR model multiplied with the `ReducInfectedSymptomsPartialImmunity` parameter.

The ratio of two variants can change over time, which affects the average transmissiblity of the disease. Infectiousness of different variants can be set in the parameters.

See *W. Koslow et al, 2022: Appropriate relaxation of non-pharmaceutical interventions minimizes the risk of a resurgence in SARS-CoV-2 infections in spite of the Delta variant* for a full description.

## Examples

The extended model is used in the paper_202107 simulation. 

Examples of the basic SECIR model can be found at:

- examples/secir.cpp
- examples/secir_ageres.cpp
- examples/parameter_study_secir.cpp
