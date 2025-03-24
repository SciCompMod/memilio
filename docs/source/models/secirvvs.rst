SECIR model with COVID-19 variants and vaccinations
=====================================================

This model extends the basic SECIR model by adding vaccinations and allowing the implicit modeling of a newly arriving variant that takes hold.

Vaccinations are modeled by adding compartments for partially and fully vaccinated persons. ``Partially and fully vaccinated`` is to be understood in this context as the person having received a first and second vaccine shot as in 2021. These model lines can be reused by resetting parameters. Persons that have recovered from the disease are treated as fully vaccinated from that time forward. Vaccinated persons are added on every day of simulation, see parameters ``DailyPartialVaccinations`` and ``DailyFullVaccinations``. All groups can get an infection or get reinfected. Vaccinated persons are less likely to develop symptoms. For example, the probability to develop symptoms when carrying the virus is the base probability from the SECIR model multiplied with the ``ReducInfectedSymptomsPartialImmunity`` parameter.

The ratio of two variants can change over time, which affects the average transmissibility of the disease. Infectiousness of different variants can be set in the parameters.

Below is an overview of the model architecture and its compartments.

.. image:: https://github.com/SciCompMod/memilio/assets/69154294/5d1b72ec-2f45-44a4-8eba-b77533c9e6cf
   :alt: SECIRVVS_model

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Mathematical variable
     - C++ variable name
     - Description
   * - :math:`\lambda_{N,i} = \rho_{N,i} \sum_j \phi_{i,j}\frac{\xi_{I_{NS}} (I_{NS,N,j} + I_{NS,PI,j} + I_{NS,II,j}) + \xi_{I_{Sy}} (I_{Sy,N,j} + I_{Sy,PI,j} + I_{Sy,II,j})}{N_j^{D^\perp}}`
     - ``ext_inf_force_dummy``
     - Force of infection for susceptibles located in the naive immunity level.
   * - :math:`\lambda_{PI,i} = \rho_{PI,i}\sum_j \phi_{i,j}\frac{\xi_{I_{NS}} (I_{NS,N,j} + I_{NS,PI,j} + I_{NS,II,j}) + \xi_{I_{Sy}} (I_{Sy,N,j} + I_{Sy,PI,j} + I_{Sy,II,j})}{N_j^{D^\perp}}`
     - ``reducExposedPartialImmunity * ext_inf_force_dummy``
     - Force of infection for susceptibles located in the partial immunity level.
   * - :math:`\lambda_{II,i} = \rho_{II}\sum_j \phi_{i,j}\frac{\xi_{I_{NS}} (I_{NS,N,j} + I_{NS,PI,j} + I_{NS,II,j}) + \xi_{I_{Sy}} (I_{Sy,N,j} + I_{Sy,PI,j} + I_{Sy,II,j})}{N_j^{D^\perp}}`
     - ``reducExposedImprovedImmunity * ext_inf_force_dummy``
     - Force of infection for susceptibles located in the improved immunity level.
   * - :math:`\phi`
     - ``ContactPatterns``
     - Matrix of daily contact rates / number of daily contacts between different age groups.
   * - :math:`\rho`
     - ``TransmissionProbabilityOnContact``
     - Transmission risk for people located in one of the susceptible compartments.
   * - :math:`\xi_{I_{NS}}`
     - ``RelativeTransmissionNoSymptoms``
     - Proportion of asymptomatically infected people who are not isolated (time-dependent if ``TestAndTraceCapacity`` is used).
   * - :math:`\xi_{I_{Sy}}`
     - ``riskFromInfectedSymptomatic``
     - Proportion of symptomatically infected people who are not isolated (time-dependent if ``TestAndTraceCapacity`` is used).
   * - :math:`N_j^{D^\perp}`
     - ``Nj``
     - Sum of all living individuals of age groups j.
   * - :math:`T_{E}`
     - ``TimeExposed``
     - Time in days an individual stays in the Exposed compartment.
   * - :math:`T_{I_{NS}}`
     - ``TimeInfectedNoSymptoms``
     - Time in days an individual stays in the InfectedNoSymptoms compartment.
   * - :math:`T_{I_{Sy}}`
     - ``TimeInfectedSymptoms``
     - Time in days an individual stays in the InfectedSymptoms compartment.
   * - :math:`T_{I_{Sev}}`
     - ``TimeInfectedSevere``
     - Time in days an individual stays in the InfectedSevere compartment.
   * - :math:`T_{I_{Cr}}`
     - ``TimeInfectedCritical``
     - Time in days an individual stays in the InfectedCritical compartment.
   * - :math:`\mu_{I_{NS}}^{I_{Sy}}`
     - ``1 - RecoveredPerInfectedNoSymptoms``
     - Probability of transition from compartment InfectedNoSymptoms to InfectedSymptoms.
   * - :math:`\mu_{I_{Sy}}^{I_{Sev}}`
     - ``SeverePerInfectedSymptoms``
     - Probability of transition from compartment InfectedSymptoms to InfectedSevere.
   * - :math:`\mu_{I_{Sev}}^{I_{Cr}}`
     - ``CriticalPerSevere``
     - Probability of transition from compartment InfectedSevere to InfectedCritical.
   * - :math:`\mu_{I_{Cr}}^{D}`
     - ``DeathsPerCritical``
     - Probability of dying when located in the InfectedCritical compartment.
   * - :math:`\kappa`
     - ``ReducTimeInfectedMild``
     - Reduction factor for time intervals for specific partial and improved immunity compartments.

Examples
--------

The extended model is used in the ``2021_vaccination_sarscov2_delta_germany`` simulation. An easier example can be found in the
`examples/ode_secirvvs.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secirvvs.cpp>`_.

Examples of the basic SECIR model can be found at:

- `examples/ode_secir.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secir.cpp>`_
- `examples/ode_secir_ageres.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secir_ageres.cpp>`_
- `examples/ode_secir_parameter_study.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secir_parameter_study.cpp>`_

Overview of the ``osecirvvs`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::osecirvvs