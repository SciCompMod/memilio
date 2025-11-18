.. include:: ../../literature.rst

ODE-based SECIRTS-type model including multi-layer waning immunity
==================================================================

This model extends the :doc:`ODE-SECIRVVS model <cpp/osecirvvs>` by adding waning immunity and introducing temporary immunity states that change the meaning of recovery.
Like the ODE-SECIRVVS model, the ODE-SECIRTS model has three layers of immunity: naive immunity, partial immunity, improved immunity.

Additionally, waning immunity is defined by the parameters ``TimeWaningPartialImmunity``, ``TimeWaningImprovedImmunity``, ``TimeTemporaryImmunityPI``, and ``TimeTemporaryImmunityII``. The parameters ``TimeWaningPartialImmunity`` and ``TimeWaningImprovedImmunity`` represent the (mean) duration after which an individual transitions from one immunity layer to the next weaker one due to waning immunity, assuming no vaccination or recovery from infection has occurred during this period. Similarly, the parameters ``TimeTemporaryImmunityPI`` and ``TimeTemporaryImmunityII`` denote the (mean) duration of temporary immunity following exposure to the virus, either through vaccination or recovery. During this state of temporary immunity, individuals are protected from reinfection and are incapable of transmitting the virus to others. Should individuals previously reside in the naive or partial immunity layer, their stay in the temporary immunity state results in a transition to the next stronger immunity layer.

For more details about the model, we refer to |Novel_travel_time_aware_metapopulation_models|.

The complete system of equations can be found in the supplementary material: `doi:10.1371/journal.pcbi.1012630.s001 <https://doi.org/10.1371/journal.pcbi.1012630.s001>`_.

Below is an overview of the model architecture and its compartments.

.. image:: https://github.com/SciCompMod/memilio/assets/69154294/6dec331f-bd91-410f-be5e-c8cf6eb0572b
   :alt: SECIRTS_model

Infection States
----------------

The model extends the ODE-SECIRVVS model by adding temporary immunity states and flow paths for waning immunity. It contains the following list of **InfectionState**\s:

.. code-block:: RST

    // Naive immunity
    `SusceptibleNaive`
    `ExposedNaive`
    `InfectedNoSymptomsNaive`
    `InfectedNoSymptomsNaiveConfirmed`
    `InfectedSymptomsNaive`
    `InfectedSymptomsNaiveConfirmed`
    `InfectedSevereNaive`
    `InfectedCriticalNaive`
    `DeadNaive`
    
    // Temporary immunity after infection in naive state
    `TemporaryImmunePartialImmunity`
    
    // Partial immunity
    `SusceptiblePartialImmunity`
    `ExposedPartialImmunity`
    `InfectedNoSymptomsPartialImmunity`
    `InfectedNoSymptomsPartialImmunityConfirmed`
    `InfectedSymptomsPartialImmunity`
    `InfectedSymptomsPartialImmunityConfirmed`
    `InfectedSeverePartialImmunity`
    `InfectedCriticalPartialImmunity`
    `DeadPartialImmunity`
    
    // Temporary immunity after infection in partial or improved immunity state
    `TemporaryImmuneImprovedImmunity`
    
    // Improved immunity
    `SusceptibleImprovedImmunity`
    `ExposedImprovedImmunity`
    `InfectedNoSymptomsImprovedImmunity`
    `InfectedNoSymptomsImprovedImmunityConfirmed`
    `InfectedSymptomsImprovedImmunity`
    `InfectedSymptomsImprovedImmunityConfirmed`
    `InfectedSevereImprovedImmunity`
    `InfectedCriticalImprovedImmunity`
    `DeadImprovedImmunity`

Infection State Transitions
---------------------------

The ODE-SECIRTS model is implemented as a **FlowModel**, which defines the derivatives of each flow between compartments. A key difference from the ODE-SECIRVVS model is that vaccinations in the ODE-SECIRTS model are implemented as flows within the ODE system rather than discrete events.

The model has the following state trnsitions:

.. code-block:: RST

    // Naive immunity flows
    `SusceptibleNaive,                            ExposedNaive`
    `SusceptibleNaive,                            TemporaryImmunePartialImmunity` (vaccination)
    `ExposedNaive,                                InfectedNoSymptomsNaive`
    `InfectedNoSymptomsNaive,                     InfectedSymptomsNaive`
    `InfectedNoSymptomsNaive,                     TemporaryImmunePartialImmunity` (recovery)
    `InfectedNoSymptomsNaiveConfirmed,            InfectedSymptomsNaiveConfirmed`
    `InfectedNoSymptomsNaiveConfirmed,            TemporaryImmunePartialImmunity` (recovery)
    `InfectedSymptomsNaive,                       InfectedSevereNaive`
    `InfectedSymptomsNaive,                       TemporaryImmunePartialImmunity` (recovery)
    `InfectedSymptomsNaiveConfirmed,              InfectedSevereNaive`
    `InfectedSymptomsNaiveConfirmed,              TemporaryImmunePartialImmunity` (recovery)
    `InfectedSevereNaive,                         InfectedCriticalNaive`
    `InfectedSevereNaive,                         TemporaryImmunePartialImmunity` (recovery)
    `InfectedSevereNaive,                         DeadNaive`
    `InfectedCriticalNaive,                       DeadNaive`
    `InfectedCriticalNaive,                       TemporaryImmunePartialImmunity` (recovery)
    
    // Temporary immunity to partial immunity transition
    `TemporaryImmunePartialImmunity,              SusceptiblePartialImmunity`
    
    // Partial immunity flows
    `SusceptiblePartialImmunity,                  ExposedPartialImmunity`
    `SusceptiblePartialImmunity,                  TemporaryImmuneImprovedImmunity` (vaccination)
    `ExposedPartialImmunity,                      InfectedNoSymptomsPartialImmunity`
    // ... similar flows for other partial immunity states leading to recovery in TemporaryImmuneImprovedImmunity ...
    
    // Temporary immunity to improved immunity transition
    `TemporaryImmuneImprovedImmunity,             SusceptibleImprovedImmunity`
    
    // Improved immunity flows
    `SusceptibleImprovedImmunity,                 ExposedImprovedImmunity`
    `SusceptibleImprovedImmunity,                 TemporaryImmuneImprovedImmunity` (booster vaccination)
    // ... similar flows for other improved immunity states leading to recovery in TemporaryImmuneImprovedImmunity ...
    
    // Waning immunity flows
    `SusceptibleImprovedImmunity,                 SusceptiblePartialImmunity`
    `SusceptiblePartialImmunity,                  SusceptibleNaive`

Sociodemographic Stratification
-------------------------------

Like the other ODE-SECIR models, the ODE-SECIRTS model can be stratified by one sociodemographic dimension, typically age groups. This stratification is important for modeling different vaccination rates, symptom severities, mortality risks, and immunity waning rates across age groups. The dimension is denoted 
**AgeGroup** but can also be used for other interpretations.
For stratifications with two or more dimensions, see :doc:`Model Creation <../ode_creation>`.

Parameters
----------

The model includes all parameters from the ODE-SECIRVVS model as well as additional parameters specific to waning and temporary immunity states:

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
   * - :math:`DV_{part,i}(t)`
     - ``DailyPartialVaccinations``
     - Daily number of first-dose vaccinations per age group.
   * - :math:`DV_{full,i}(t)`
     - ``DailyFullVaccinations``
     - Daily number of second-dose vaccinations per age group.
   * - :math:`DV_{boost,i}(t)`
     - ``DailyBoosterVaccinations``
     - Daily number of booster vaccinations per age group.
   * - :math:`T_{V1}`
     - ``DaysUntilEffectivePartialVaccination``
     - Time in days until first vaccine dose takes full effect.
   * - :math:`T_{V2}`
     - ``DaysUntilEffectiveImprovedVaccination``
     - Time in days until second vaccine dose takes full effect.
   * - :math:`T_{V3}`
     - ``DaysUntilEffectiveBoosterImmunity``
     - Time in days until booster vaccine dose takes full effect.
   * - :math:`s`
     - ``Seasonality``
     - Seasonal variation factor affecting transmission.
   * - :math:`ICU_{cap}`
     - ``ICUCapacity``
     - Maximum ICU capacity in the modeled region.
   * - :math:`\nu`
     - ``InfectiousnessNewVariant``
     - Relative infectiousness of a new variant compared to the original strain.
   * - :math:`t_{newVar}`
     - ``StartDayNewVariant``
     - Day on which the new variant is introduced in the simulation.

Initial conditions
------------------

The initial conditions of the model are represented by the class **Populations** which defines the number of individuals in each sociodemographic group and **InfectionState**. Before running a simulation, the initial values for each compartment across all immunity levels have to be set. This can be done via:

.. code-block:: cpp

    for (mio::AgeGroup i = 0; i < nb_groups; i++) {
        // population
        model.populations[{i, mio::osecirts::InfectionState::ExposedNaive}]                                = 20;
        model.populations[{i, mio::osecirts::InfectionState::ExposedImprovedImmunity}]                     = 20;
        model.populations[{i, mio::osecirts::InfectionState::ExposedPartialImmunity}]                      = 20;
        model.populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsNaive}]                     = 30;
        model.populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsNaiveConfirmed}]            = 0;
        model.populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunity}]           = 30;
        model.populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]  = 0;
        model.populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunity}]          = 30;
        model.populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 0;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSymptomsNaive}]                       = 40;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSymptomsNaiveConfirmed}]              = 0;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity}]             = 40;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]    = 0;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunity}]            = 40;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]   = 0;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSevereNaive}]                         = 30;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSevereImprovedImmunity}]              = 30;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSeverePartialImmunity}]               = 30;
        model.populations[{i, mio::osecirts::InfectionState::InfectedCriticalNaive}]                       = 20;
        model.populations[{i, mio::osecirts::InfectionState::InfectedCriticalPartialImmunity}]             = 20;
        model.populations[{i, mio::osecirts::InfectionState::InfectedCriticalImprovedImmunity}]            = 20;
        model.populations[{i, mio::osecirts::InfectionState::SusceptibleNaive}]                            = 1000;
        model.populations[{i, mio::osecirts::InfectionState::SusceptiblePartialImmunity}]                  = 1200;
        model.populations[{i, mio::osecirts::InfectionState::SusceptibleImprovedImmunity}]                 = 1000;
        model.populations[{i, mio::osecirts::InfectionState::TemporaryImmunePartialImmunity}]              = 60;
        model.populations[{i, mio::osecirts::InfectionState::TemporaryImmuneImprovedImmunity}]             = 70;
        model.populations[{i, mio::osecirts::InfectionState::DeadNaive}]                                   = 0;
        model.populations[{i, mio::osecirts::InfectionState::DeadPartialImmunity}]                         = 0;
        model.populations[{i, mio::osecirts::InfectionState::DeadImprovedImmunity}]                        = 0;
    }

After setting the initial populations, the daily vaccination parameters, which are directly integrated into the ODE system in this model, also need to be configured:

.. code-block:: cpp

    const size_t daily_vaccinations = 10;
    const size_t num_days = 300;
    model.parameters.get<mio::osecirts::DailyPartialVaccinations<double>>().resize(mio::SimulationDay(num_days));
    model.parameters.get<mio::osecirts::DailyFullVaccinations<double>>().resize(mio::SimulationDay(num_days));
    model.parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>().resize(mio::SimulationDay(num_days));
    for (size_t i = 0; i < num_days; ++i) {
        for (mio::AgeGroup j = 0; j < nb_groups; ++j) {
            auto num_vaccinations = static_cast<double>(i * daily_vaccinations);
            model.parameters.get<mio::osecirts::DailyPartialVaccinations<double>>()[{j, mio::SimulationDay(i)}] =
                num_vaccinations;
            model.parameters.get<mio::osecirts::DailyFullVaccinations<double>>()[{j, mio::SimulationDay(i)}] =
                num_vaccinations;
            model.parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>()[{j, mio::SimulationDay(i)}] =
                num_vaccinations;
        }
    }

Nonpharmaceutical Interventions
-------------------------------

Like the other ODE-SECIR models, the ODE-SECIRTS model supports nonpharmaceutical interventions (NPIs) through dampings in the contact matrix. These dampings reduce the contact rates between different groups to simulate interventions like lockdowns.

Basic dampings can be added to the contact matrix as follows:

.. code-block:: cpp

    // Create a contact matrix with baseline contact rates
    auto& contacts = model.parameters.get<mio::osecirts::ContactPatterns<double>>();
    auto& contact_matrix = contacts.get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(0.5);
    contact_matrix[0].get_baseline().diagonal().setConstant(5.0);
    
    // Add a damping that reduces contacts by 30% starting at day 5
    contact_matrix[0].add_damping(0.3, mio::SimulationTime(5.0));

The model also supports dynamic NPIs based on epidemic thresholds:

.. code-block:: cpp
    
    // Set threshold-based triggers for NPIs
    auto& dynamic_npis = model.parameters.get<mio::osecirts::DynamicNPIsInfectedSymptoms<double>>();
    dynamic_npis.set_interval(mio::SimulationTime(3.0));  // Check every 3 days
    dynamic_npis.set_duration(mio::SimulationTime(14.0)); // Apply for 14 days
    dynamic_npis.set_base_value(100'000);                // Per 100,000 population
    dynamic_npis.set_threshold(200.0, dampings);         // Trigger at 200 cases per 100,000

For more complex scenarios, such as real-world lockdown modeling, detailed NPIs with location-specific dampings can be implemented. For further details, see the documentation of the ODE-SECIR model. 


Simulation
----------

The ODE-SECIRTS model offers the same simulation functions as the other ODE-SECIR models:

1. **simulate**: Standard simulation that tracks the compartment sizes over time
2. **simulate_flows**: Extended simulation that additionally tracks the flows between compartments

Standard simulation:

.. code-block:: cpp

    double t0 = 0;       // Start time
    double tmax = 50;    // End time
    double dt = 0.1;     // Time step
    
    // Run a standard simulation
    mio::TimeSeries<double> result = mio::osecirts::simulate<double>(t0, tmax, dt, model);

During simulation, the model handles several special processes:

1. **Vaccinations**: Unlike the ODE-SECIRVVS model, vaccinations are integrated directly into the ODE system through flows from susceptible compartments to temporary immunity compartments.

2. **Variant Evolution**: The `apply_variant` function updates the transmission probability based on the existance of a new variant over time, similar to other ODE-SECIR models.

For both simulation types, you can also specify a custom integrator:

.. code-block:: cpp

    auto integrator = std::make_unique<mio::RKIntegratorCore>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    
    mio::TimeSeries<double> result = mio::osecirts::simulate(t0, tmax, dt, model, std::move(integrator));

Output
------

The output of the simulation is a `mio::TimeSeries` object containing the sizes of each compartment at each time point. For a standard simulation, you can access the results as follows:

.. code-block:: cpp

    // Get the number of time points
    auto num_points = static_cast<size_t>(result.get_num_time_points());
    
    // Access data at a specific time point
    Eigen::VectorXd value_at_time_i = result.get_value(i);
    double time_i = result.get_time(i);
    
    // Access the last time point
    Eigen::VectorXd last_value = result.get_last_value();

You can print the simulation results as a formatted table:

.. code-block:: cpp

    // Print results to console with default formatting
    result.print_table();

The order of the compartments follows the definition in the `InfectionState` enum.

Additionally, you can export the results to a CSV file for further analysis or visualization:

.. code-block:: cpp

    // Export results to CSV
    result.export_csv("simulation_results.csv");

Visualization
-------------

To visualize the results of a simulation, you can use the Python package :doc:`m-plot <../../python/m-plot>`
and its documentation.

Examples
--------

To get started with the ODE-SECIRTS model, check out the code example in the MEmilio repository:
`examples/ode_secirts.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secirts.cpp>`_.

Overview of the ``osecirts`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::osecirts
