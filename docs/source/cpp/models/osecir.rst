ODE-based SECIR-type model
===========================

The ODE-based SECIR-type module models and simulates the epidemic using an ODE-based SECIR-type model approach. 
The model is particularly suited for pathogen with pre- or asymptomatic infection states and when severe or critical
states are possible. The model assumes perfect immunity after recovery and is thus only suited for epidemic use cases.
In the following, we present the model in detail.

The infection states and the transitions (also see next to sections) are visualized in the following image.

.. image:: https://github.com/SciCompMod/memilio/assets/70579874/46b09e8a-d083-4ef9-8328-21975890b60f
   :alt: secir_model

Infection States
----------------

The model contains the following list of **InfectionState**\s:

.. code-block:: RST

    `Susceptible`
    `Exposed`
    `InfectedNoSymptoms`
    `InfectedNoSymptomsConfirmed`
    `InfectedSymptoms`
    `InfectedSymptomsConfirmed`
    `InfectedSevere`
    `InfectedCritical`
    `Recovered`
    `Dead`

While the states `InfectedNoSymptomsConfirmed` and `InfectedSymptomsConfirmed` are available, they are not used in the 
current implementation and detection is only modeled implicitly through detection rates based on user-defined criteria.

Infection State Transitions
---------------------------

The ODE-SECIR model is implemented as **FlowModel**. With just minimal overhead, the **FlowModel** computes the new 
transmissions, infections, and hospitalizations explicitly in every time step instead of only computing the aggregated 
compartment values. The defined transitions `FromState, ToState` are:

.. code-block:: RST
    `Susceptible,                 Exposed`
    `Exposed,                     InfectedNoSymptoms`
    `InfectedNoSymptoms,          InfectedSymptoms`
    `InfectedNoSymptoms,          Recovered`
    `InfectedNoSymptomsConfirmed, InfectedSymptomsConfirmed`
    `InfectedNoSymptomsConfirmed, Recovered`
    `InfectedSymptoms,            InfectedSevere`
    `InfectedSymptoms,            Recovered`
    `InfectedSymptomsConfirmed,   InfectedSevere`
    `InfectedSymptomsConfirmed,   Recovered`
    `InfectedSevere,              InfectedCritical`
    `InfectedSevere,              Recovered`
    `InfectedSevere,              Dead`
    `InfectedCritical,            Dead`
    `InfectedCritical,            Recovered`


Sociodemographic Stratification
-------------------------------

In the ODE-SECIR model, the population can be stratified by one sociodemographic dimension. This dimension is denoted 
**AgeGroup** but can also be used for other interpretations. For stratifications with two or more dimensions, 
see :doc:`Model Creation <../ode_creation>`.


Parameters
----------

The model implements the following parameters.

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Mathematical variable
     - C++ variable name
     - Description
   * - :math:`\phi`
     - ``ContactPatterns``
     - Matrix of daily contact rates / number of daily contacts between different age groups.
   * - :math:`\rho`
     - ``TransmissionProbabilityOnContact``
     - Transmission risk for people located in one of the susceptible compartments.
   * - :math:`\xi_{I_{NS}}`
     - ``RelativeTransmissionNoSymptoms``
     - Proportion of nonsymptomatically infected people who are not isolated.
   * - :math:`\xi_{I_{Sy}}`
     - ``riskFromInfectedSymptomatic``
     - Proportion of infected people with symptoms who are not isolated (time-dependent if ``TestAndTraceCapacity`` used).
   * - :math:`N_j`
     - ``Nj``
     - Total population of age group :math:`j`.
   * - :math:`D_i`
     - ``Di``
     - Number of deaths in age group :math:`i`.
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
     - Probability of dying when in compartment InfectedCritical.


Initial conditions
------------------

TODO...


.. _Nonpharmaceutical Interventions:
Nonpharmaceutical Interventions
-------------------------------

TODO...


Simulation
----------

TODO...


Output
------

TODO...


Visualization
-------------

TODO...

    
Examples
--------

Different examples can be found at:

- `examples/ode_secir.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secir.cpp>`_
- `examples/ode_secir_ageres.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secir_ageres.cpp>`_
- `examples/ode_secir_parameter_study.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secir_parameter_study.cpp>`_  

Overview of the ``osecir`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::osecir




OLD Structure
---------


3. **Dampings**: A ``Damping`` object is the combination of a particular day and a multiplicative factor that changes the contact patterns. Dampings can be overwritten by or combined with dampings at later times. In order to avoid discontinuities in pattern changes, the transition is smoothed by a cosine S-type function over an interval of maximum length one day. The resulting contact rate satisfies C^1-smoothness between two consecutive dampings.
4. **SECIR**: Implements an *age-resolved ODE-model*, based on the non-age-resolved model as described in
   `https://www.medrxiv.org/content/10.1101/2020.04.04.20053637v1 <https://www.medrxiv.org/content/10.1101/2020.04.04.20053637v1>`_. It uses the compartments ``Susceptible (S)``, ``Exposed (E)``, ``InfectedNoSymptoms (I_NS)``, ``InfectedSymptoms (I_Sy)``, ``InfectedSevere (I_Sev)``, ``InfectedCritical (I_Cr)``, ``Recovered (R)`` and ``Dead``. The final model has been published in
   `https://doi.org/10.1016/j.mbs.2021.108648 <https://doi.org/10.1016/j.mbs.2021.108648>`_. ``Exposed`` individuals are not infectious and it is assumed that subsequent infections can occur at random during the infectious period before the onset of symptoms, i.e., from state ``InfectedNoSymptoms (I_NS)``. Recovered people remain immune. Severely or critically infected individuals are assumed to be isolated. The ``Model`` uses ``Populations`` to model different 'groups' of a particular age range (first dimension) and an ``InfectionState`` (second dimension). Parameters are set as ``Parameters``; they contain contact patterns in the form of an ``UncertainContactMatrix`` and an extended set of pathogen-dependent parameters.
5. **Parameter Space**: Factory class for the ``Parameters`` to set distributions for the different parameters and to provide the opportunity to sample from this parameter set containing random distributions.
6. **Parameter Studies**: Method to be called on a set of ``Parameters`` with a given set of random distributions to sample from the distributions and run ensemble simulations with the obtained samples.


Simulation
----------

The simulation runs in discrete time steps using a numerical integration scheme. At each time step, a portion of the population in each age-aware compartment moves from the current compartment to a new one. Different numerical integration schemes are available (see the ``math`` folder). The Simulation class handles the parameters and the numerical integrator, and it also stores the result. Ensemble runs can be performed using the Parameter Studies class as soon as random distributions are set for all parameters. This can be done using the Parameter Space class.




