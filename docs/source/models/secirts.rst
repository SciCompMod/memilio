SECIRTS-type model including multi-layer waning immunity
=========================================================

This model extends the SECIRVVS model by adding waning immunity and introducing temporary immunity states that change the meaning of recovery.
We structured the model into three layers of immunity: naive immunity, partial immunity, improved immunity.

In the model, waning immunity is defined by the parameters ``TimeWaningPartialImmunity``, ``TimeWaningImprovedImmunity``, ``TimeTemporaryImmunityPI``, and ``TimeTemporaryImmunityII``. The parameters ``TimeWaningPartialImmunity`` and ``TimeWaningImprovedImmunity`` represent the (mean) duration after which an individual transitions from one immunity layer to the next less protected one due to waning immunity, assuming no vaccination or recovery from infection has occurred during this period. Also, the parameters ``TimeTemporaryImmunityPI`` and ``TimeTemporaryImmunityII`` denote the (mean) duration of temporary immunity following exposure to the virus, either through vaccination or recovery. During this state of temporary immunity, individuals are protected from reinfection and are incapable of transmitting the virus to others. Should individuals previously reside in the naive or partial immunity layer, their stay in the temporary immunity state results in a transition to the next more protected immunity layer.

For more details about the model, we refer to `1 <https://www.medrxiv.org/content/10.1101/2024.03.01.24303602v3>`_.

Below is an overview of the model architecture and its compartments.

.. image:: https://github.com/SciCompMod/memilio/assets/69154294/6dec331f-bd91-410f-be5e-c8cf6eb0572b
   :alt: SECIRTS_model

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Mathematical variable
     - C++ variable name
     - Description
   * - :math:`\lambda_{N,i} = \rho_{N,i} \sum_j \phi_{i,j}\frac{\xi_{I_{NS}} \Bigl(I_{NS,N,j} + I_{NS,PI,j} + I_{NS,II,j}\Bigr) + \xi_{I_{Sy}} \Bigl(I_{Sy,N,j} + I_{Sy,PI,j} + I_{Sy,II,j}\Bigr)}{N_j^{D^\perp}}`
     - ``ext_inf_force_dummy``
     - Force of infection for susceptibles located in the naive immunity level.
   * - :math:`\lambda_{PI,i} = \rho_{PI,i}\sum_j \phi_{i,j}\frac{\xi_{I_{NS}} \Bigl(I_{NS,N,j} + I_{NS,PI,j} + I_{NS,II,j}\Bigr) + \xi_{I_{Sy}} \Bigl(I_{Sy,N,j} + I_{Sy,PI,j} + I_{Sy,II,j}\Bigr)}{N_j^{D^\perp}}`
     - ``reducExposedPartialImmunity * ext_inf_force_dummy``
     - Force of infection for susceptibles located in the partial immunity level.
   * - :math:`\lambda_{II,i} = \rho_{II}\sum_j \phi_{i,j}\frac{\xi_{I_{NS}} \Bigl(I_{NS,N,j} + I_{NS,PI,j} + I_{NS,II,j}\Bigr) + \xi_{I_{Sy}} \Bigl(I_{Sy,N,j} + I_{Sy,PI,j} + I_{Sy,II,j}\Bigr)}{N_j^{D^\perp}}`
     - ``reducExposedImprovedImmunity * ext_inf_force_dummy``
     - Force of infection for susceptibles located in the improved immunity level.
   * - :math:`\phi`
     - ``ContactPatterns``
     - Matrix of daily contact rates, i.e., number of daily contacts between different age groups.
   * - :math:`\rho`
     - ``TransmissionProbabilityOnContact``
     - Transmission risk for people located in one of the susceptible compartments.
   * - :math:`\xi_{I_{NS}}`
     - ``RelativeTransmissionNoSymptoms``
     - Proportion of asymptomatically infected people who are not isolated (time-dependent if ``TestAndTraceCapacity`` used).
   * - :math:`\xi_{I_{Sy}}`
     - ``riskFromInfectedSymptomatic``
     - Proportion of symptomatically infected people who are not isolated (time-dependent if ``TestAndTraceCapacity`` used).
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
   * - :math:`T_{\mathcal{I}_{PI}}`
     - ``TimeTemporaryImmunityPI``
     - Time in days an individual stays in the TemporaryImmunePartialImmunity compartment.
   * - :math:`T_{\mathcal{I}_{PI}}`
     - ``TimeTemporaryImmunityII``
     - Time in days an individual stays in the TemporaryImmuneImprovedImmunity compartment.
   * - :math:`T_{W_{PI}}`
     - ``TimeWaningPartialImmunity``
     - Time in days an individual stays in the SusceptiblePartialImmunity compartment before waning to the SusceptibleNaive compartment assuming no exposure occurred during this period.
   * - :math:`T_{W_{II}}`
     - ``TimeWaningImprovedImmunity``
     - Time in days an individual stays in the SusceptibleImprovedImmunity compartment before waning to the SusceptiblePartialImmunity compartment assuming no exposure occurred during this period.
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
     - Probability of dying when located in compartment InfectedCritical.
   * - :math:`\kappa`
     - ``ReducTimeInfectedMild``
     - Reduction factor for time intervals for specific partial and improved immunity compartments.


Overview of the ``osecirts`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::osecirts