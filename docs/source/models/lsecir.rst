LCT SECIR model
===============

This model is based on the Linear Chain Trick.

The Linear Chain Trick provides the option to use Erlang-distributed stay times in the compartments through the use of subcompartments. The normal ODE models have (possibly unrealistic) exponentially distributed stay times. The LCT model can still be described by an ordinary differential equation system.

For the concept see:

- `Lena Plötzke, "Der Linear Chain Trick in der epidemiologischen Modellierung als Kompromiss zwischen gewöhnlichen und Integro-Differentialgleichungen", 2023 <https://elib.dlr.de/200381/>`_ (German only)
- `P. J. Hurtado und A. S. Kirosingh, "Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models“, 2019 <https://doi.org/10.1007/s00285-019-01412-w>`_

The eight compartments

- **Susceptible** (:math:`S`), may become exposed at any time.
- **Exposed** (:math:`E`), becomes infected after some time.
- **InfectedNoSymptoms** (:math:`I_{NS}`), becomes InfectedSymptoms or Recovered after some time.
- **InfectedSymptoms** (:math:`I_{Sy}`), becomes InfectedSevere or Recovered after some time.
- **InfectedSevere** (:math:`I_{Sev}`), becomes InfectedCritical or Recovered after some time.
- **InfectedCritical** (:math:`I_{Cr}`), becomes Recovered or Dead after some time.
- **Recovered** (:math:`R`)
- **Dead** (:math:`D`)

are used to simulate the spread of the disease. It is possible to include subcompartments for the five compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere, and InfectedCritical. You can divide the population according to different groups (e.g. AgeGroups or gender) and choose parameters accordingly.

Below is an overview of the model architecture and its compartments without a stratification according to groups.

.. image:: https://github.com/SciCompMod/memilio/assets/70579874/6a5d5a95-20f9-4176-8894-c091bd48bfb7
   :alt: tikzLCTSECIR

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Mathematical variable
     - C++ variable name
     - Description
   * - :math:`\phi`
     - ``ContactPatterns``
     - Average number of contacts of a person per day.
   * - :math:`\rho`
     - ``TransmissionProbabilityOnContact``
     - Transmission risk for people located in the susceptible compartments.
   * - :math:`\xi_{I_{NS}}`
     - ``RelativeTransmissionNoSymptoms``
     - Proportion of nonsymptomatically infected people who are not isolated.
   * - :math:`\xi_{I_{Sy}}`
     - ``RiskOfInfectionFromSymptomatic``
     - Proportion of infected people with symptoms who are not isolated.
   * - :math:`N`
     - ``m_N0``
     - Total population.
   * - :math:`D`
     - ``D``
     - Number of death people.
   * - :math:`n_E`
     - Defined in ``LctStates``
     - Number of subcompartments of the Exposed compartment.
   * - :math:`n_{NS}`
     - Defined in ``LctStates``
     - Number of subcompartments of the InfectedNoSymptoms compartment.
   * - :math:`n_{Sy}`
     - Defined in ``LctStates``
     - Number of subcompartments of the InfectedSymptoms compartment.
   * - :math:`n_{Sev}`
     - Defined in ``LctStates``
     - Number of subcompartments of the InfectedSevere compartment.
   * - :math:`n_{Cr}`
     - Defined in ``LctStates``
     - Number of subcompartments of the InfectedCritical compartment.
   * - :math:`T_E`
     - ``TimeExposed``
     - Average time in days an individual stays in the Exposed compartment.
   * - :math:`T_{I_{NS}}`
     - ``TimeInfectedNoSymptoms``
     - Average time in days an individual stays in the InfectedNoSymptoms compartment.
   * - :math:`T_{I_{Sy}}`
     - ``TimeInfectedSymptoms``
     - Average time in days an individual stays in the InfectedSymptoms compartment.
   * - :math:`T_{I_{Sev}}`
     - ``TimeInfectedSevere``
     - Average time in days an individual stays in the InfectedSevere compartment.
   * - :math:`T_{I_{Cr}}`
     - ``TimeInfectedCritical``
     - Average time in days an individual stays in the InfectedCritical compartment.
   * - :math:`\mu_{I_{NS}}^{R}`
     - ``RecoveredPerInfectedNoSymptoms``
     - Probability of transition from compartment InfectedNoSymptoms to Recovered.
   * - :math:`\mu_{I_{Sy}}^{I_{Sev}}`
     - ``SeverePerInfectedSymptoms``
     - Probability of transition from compartment InfectedSymptoms to InfectedSevere.
   * - :math:`\mu_{I_{Sev}}^{I_{Cr}}`
     - ``CriticalPerSevere``
     - Probability of transition from compartment InfectedSevere to InfectedCritical.
   * - :math:`\mu_{I_{Cr}}^{D}`
     - ``DeathsPerCritical``
     - Probability of dying when in compartment InfectedCritical.

The notation of the compartments with indices here stands for subcompartments and not for age groups. Accordingly, :math:`I_{NS,n_{NS}}`, for example, stands for the number of people in the :math:`n_{NS}`-th subcompartment of the InfectedNoSymptoms compartment.

Examples
--------

A simple example can be found at the
`LCT minimal example <https://github.com/SciCompMod/memilio/tree/main/cpp/examples/lct_secir.cpp>`_.

Initialization
--------------

- The file `parameters_io <https://github.com/SciCompMod/memilio/blob/main/cpp/models/lct_secir/parameters_io.h>`_ provides functionality to compute an initial value vector for the LCT-SECIR model based on real data.

- The file `initializer_flows <https://github.com/SciCompMod/memilio/blob/main/cpp/models/lct_secir/initializer_flows.h>`_ provides functionality to compute an initial value vector for the LCT-SECIR model based on initial data in the form of a TimeSeries of InfectionTransitions. For the concept of the InfectionTransitions or flows, see also the IDE-SECIR model. This method can be particularly useful if a comparison is to be made with an IDE model with matching initialization or if the real data is in the form of flows.

Overview of the ``lsecir`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::lsecir