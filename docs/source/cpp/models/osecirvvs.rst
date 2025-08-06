SECIR model with COVID-19 variants and vaccinations
=====================================================

This model extends the basic SECIR model by adding vaccinations and allowing the implicit modeling of a newly arriving variant that takes hold.

Vaccinations are modeled by adding compartments for partially and fully vaccinated persons. ``Partially and fully vaccinated`` is to be understood in this context as the person having received a first and second vaccine shot as in 2021. These model lines can be reused by resetting parameters. Persons that have recovered from the disease are treated as fully vaccinated from that time forward. Vaccinated persons are added on every day of simulation, see parameters ``DailyPartialVaccinations`` and ``DailyFullVaccinations``. All groups can get an infection or get reinfected. Vaccinated persons are less likely to develop symptoms. For example, the probability to develop symptoms when carrying the virus is the base probability from the SECIR model multiplied with the ``ReducInfectedSymptomsPartialImmunity`` parameter.

The ratio of two variants can change over time, which affects the average transmissibility of the disease. Infectiousness of different variants can be set in the parameters.

Below is an overview of the model architecture and its compartments.

.. image:: https://github.com/SciCompMod/memilio/assets/69154294/5d1b72ec-2f45-44a4-8eba-b77533c9e6cf
   :alt: SECIRVVS_model

Infection States
----------------

The model extends the basic SECIR model by dividing the compartments based on immunity levels. It contains the following list of **InfectionState**\s:

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

All compartments with the same base state (e.g., ExposedNaive, ExposedPartialImmunity, ExposedImprovedImmunity) represent the same disease stage but with different immunity levels, which affects the later course of infection.

Infection State Transitions
---------------------------

The ODE-SECIRVVS model is implemented as a **FlowModel**, which computes the flows between compartments explicitly. The model follows the same flow pattern as the basic SECIR model but with three parallel sets of compartments representing different immunity levels.

The key characteristic of this model is that recovered individuals always end up in the improved immunity level, regardless of their starting immunity level. This represents the immunity gained after infection.

For each immunity level (Naive, PartialImmunity, ImprovedImmunity), the following transitions are defined:

.. code-block:: RST
  
    `Susceptible*,                 Exposed*`
    `Exposed*,                     InfectedNoSymptoms*`
    `InfectedNoSymptoms*,          InfectedSymptoms*`
    `InfectedNoSymptoms*,          SusceptibleImprovedImmunity` (recovery always leads to improved immunity)
    `InfectedNoSymptomsConfirmed*, InfectedSymptomsConfirmed*`
    `InfectedNoSymptomsConfirmed*, SusceptibleImprovedImmunity`
    `InfectedSymptoms*,            InfectedSevere*`
    `InfectedSymptoms*,            SusceptibleImprovedImmunity`
    `InfectedSymptomsConfirmed*,   InfectedSevere*`
    `InfectedSymptomsConfirmed*,   SusceptibleImprovedImmunity`
    `InfectedSevere*,              InfectedCritical*`
    `InfectedSevere*,              SusceptibleImprovedImmunity`
    `InfectedSevere*,              Dead*`
    `InfectedCritical*,            Dead*`
    `InfectedCritical*,            SusceptibleImprovedImmunity`

Where * stands for the immunity level suffix (Naive, PartialImmunity, or ImprovedImmunity).

**Important:** Vaccinations are not implemented as flows between compartments but are handled discretely by the simulation. At the beginning of each simulated day, susceptible individuals are moved between immunity levels according to the specified daily vaccination parameters. This discrete process is separate from ODE system and is managed by the `apply_vaccination` function in the model specific Simulation class.

Sociodemographic Stratification
-------------------------------

Like the basic SECIR model, the SECIRVVS model can be stratified by one sociodemographic dimension, typically age groups. This stratification is important for modeling different vaccination rates, symptom severities, and mortality risks across age groups.
For stratifications with two or more dimensions, see :doc:`Model Creation <../ode_creation>`.

Parameters
----------

The model includes all parameters from the basic SECIR model plus additional parameters specific to vaccination and variant modeling. Here is a comprehensive list of the key parameters:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Mathematical variable
     - C++ variable name
     - Description
   * - :math:`\phi`
     - ``ContactPatterns``
     - Matrix of daily contact rates between different age groups.
   * - :math:`\rho`
     - ``TransmissionProbabilityOnContact``
     - Transmission risk when susceptible and infected individuals make contact.
   * - :math:`\xi_{I_{NS}}`
     - ``RelativeTransmissionNoSymptoms``
     - Proportion of non-symptomatically infected people who are not isolated.
   * - :math:`\xi_{I_{Sy}}`
     - ``RiskOfInfectionFromSymptomatic``
     - Proportion of infected people with symptoms who are not isolated.
   * - :math:`\mu_{I_{NS}}^{I_{Sy}}`
     - ``1 - RecoveredPerInfectedNoSymptoms``
     - Probability of transition from InfectedNoSymptoms to InfectedSymptoms.
   * - :math:`\mu_{I_{Sy}}^{I_{Sev}}`
     - ``SeverePerInfectedSymptoms``
     - Probability of transition from InfectedSymptoms to InfectedSevere.
   * - :math:`\mu_{I_{Sev}}^{I_{Cr}}`
     - ``CriticalPerSevere``
     - Probability of transition from InfectedSevere to InfectedCritical.
   * - :math:`\mu_{I_{Cr}}^{D}`
     - ``DeathsPerCritical``
     - Probability of dying when in InfectedCritical.
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
   * - :math:`s`
     - ``Seasonality``
     - Seasonal variation factor affecting transmission.
   * - :math:`ICU_{cap}`
     - ``ICUCapacity``
     - Maximum ICU capacity in the modeled region.
   * - :math:`TTC_{cap}`
     - ``TestAndTraceCapacity``
     - Capacity to test and trace contacts of infected per day.
   * - :math:`TTC_{maxNoSym}`
     - ``TestAndTraceCapacityMaxRiskNoSymptoms``
     - Multiplier for test and trace capacity for cases without symptoms.
   * - :math:`TTC_{maxSym}`
     - ``TestAndTraceCapacityMaxRiskSymptoms``
     - Multiplier for test and trace capacity for symptomatic cases.
   * - :math:`T_{dyndelay}`
     - ``DynamicNPIsImplementationDelay``
     - Delay in days for implementing dynamic NPIs after threshold exceedance.
   * - :math:`\lambda_{N,i}`
     - ``ext_inf_force_dummy``
     - Force of infection for susceptibles with naive immunity.
   * - :math:`\lambda_{PI,i}`
     - ``reducExposedPartialImmunity * ext_inf_force_dummy``
     - Force of infection for susceptibles with partial immunity.
   * - :math:`\lambda_{II,i}`
     - ``reducExposedImprovedImmunity * ext_inf_force_dummy``
     - Force of infection for susceptibles with improved immunity.
   * - :math:`DV_{PI,i}`
     - ``DailyPartialVaccinations``
     - Daily number of first-dose vaccinations per age group.
   * - :math:`DV_{II,i}`
     - ``DailyFullVaccinations``
     - Daily number of second-dose vaccinations per age group.
   * - :math:`T_{vaccGap}`
     - ``VaccinationGap``
     - Time in days between first and second vaccine dose.
   * - :math:`T_{V1}`
     - ``DaysUntilEffectivePartialImmunity``
     - Time in days until first vaccine dose takes full effect.
   * - :math:`T_{V2}`
     - ``DaysUntilEffectiveImprovedImmunity``
     - Time in days until second vaccine dose takes full effect.
   * - :math:`\delta_{E,PI}`
     - ``ReducExposedPartialImmunity``
     - Factor to reduce infection risk for persons with partial immunity.
   * - :math:`\delta_{E,II}`
     - ``ReducExposedImprovedImmunity``
     - Factor to reduce infection risk for persons with improved immunity.
   * - :math:`\delta_{Sy,PI}`
     - ``ReducInfectedSymptomsPartialImmunity``
     - Factor to reduce risk of developing symptoms for persons with partial immunity.
   * - :math:`\delta_{Sy,II}`
     - ``ReducInfectedSymptomsImprovedImmunity``
     - Factor to reduce risk of developing symptoms for persons with improved immunity.
   * - :math:`\delta_{Sev,PI}`
     - ``ReducInfectedSevereCriticalDeadPartialImmunity``
     - Factor to reduce hospitalization/mortality risk for partial immunity.
   * - :math:`\delta_{Sev,II}`
     - ``ReducInfectedSevereCriticalDeadImprovedImmunity``
     - Factor to reduce hospitalization/mortality risk for improved immunity.
   * - :math:`\kappa`
     - ``ReducTimeInfectedMild``
     - Reduction factor for infectious periods with partial or improved immunity.
   * - :math:`\nu`
     - ``InfectiousnessNewVariant``
     - Relative infectiousness of a new variant compared to the original variant.
   * - :math:`t_{newVar}`
     - ``StartDayNewVariant``
     - Day on which the new variant is introduced in the simulation.

Initial conditions
------------------

The initial conditions of the model are represented by the class **Populations** which defines the number of individuals in each sociodemographic group and **InfectionState**. Before running a simulation, you should set the initial values for each compartment across all immunity levels.

Below is a example showing how to initialize all compartments for the SECIRVVS model:

.. code-block:: cpp

    for (mio::AgeGroup i = 0; i < model.parameters.get_num_groups(); i++) {
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedNaive}]                                = 10;
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedImprovedImmunity}]                     = 11;
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedPartialImmunity}]                      = 12;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}]                     = 13;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}]            = 13;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}]           = 14;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]  = 14;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}]          = 15;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 15;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}]                       = 5;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}]              = 5;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}]             = 6;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]    = 6;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}]            = 7;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]   = 7;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereNaive}]                         = 8;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereImprovedImmunity}]              = 1;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSeverePartialImmunity}]               = 2;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalNaive}]                       = 3;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalPartialImmunity}]             = 4;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalImprovedImmunity}]            = 5;
        model.populations[{i, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}]                 = 6;
        model.populations[{i, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}]                  = 7;
        model.populations[{i, mio::osecirvvs::InfectionState::DeadNaive}]                                   = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::DeadPartialImmunity}]                         = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::DeadImprovedImmunity}]                        = 0;
        
        // Set the SusceptibleNaive compartment as the difference from the total population
        model.populations.set_difference_from_group_total<mio::AgeGroup>(
            {i, mio::osecirvvs::InfectionState::SusceptibleNaive}, 1000);
    }

After setting the initial populations, you also need to set the vaccination parameters. Vaccinations are realized discretely by moving susceptible individuals to the partially and fully vaccinated immunity levels at the beginning of each day of the simulation. For example, to set daily vaccinations that increase over time:

.. code-block:: cpp
    
    // Prepare and resize vaccinations parameter for the entire simulation period
    const size_t daily_vaccinations = 10;
    model.parameters.get<mio::osecirvvs::DailyPartialVaccinations<double>>().resize(
        mio::SimulationDay((size_t)tmax + 1));
    model.parameters.get<mio::osecirvvs::DailyFullVaccinations<double>>().resize(
        mio::SimulationDay((size_t)tmax + 1));
        
    // Set increasing number of vaccination over time
    for (size_t i = 0; i < tmax + 1; ++i) {
        auto num_vaccinations = static_cast<double>(i * daily_vaccinations);
        model.parameters
            .get<mio::osecirvvs::DailyPartialVaccinations<double>>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
        model.parameters
            .get<mio::osecirvvs::DailyFullVaccinations<double>>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
    }

.. _Nonpharmaceutical Interventions:
Nonpharmaceutical Interventions
-------------------------------

The SECIRVVS model supports nonpharmaceutical interventions (NPIs) through dampings in the contact matrix. These dampings reduce the contact rates between different groups to simulate interventions like lockdowns.

Basic dampings can be added to the contact matrix as follows:

.. code-block:: cpp

    // Create a contact matrix with baseline contact rates
    auto& contacts = model.parameters.get<mio::osecirvvs::ContactPatterns<double>>();
    auto& contact_matrix = contacts.get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(0.5);
    contact_matrix[0].get_baseline().diagonal().setConstant(5.0);
    
    // Add a damping that reduces contacts by 30% starting at day 5
    contact_matrix[0].add_damping(0.3, mio::SimulationTime(5.0));

For more complex scenarios, such as real-world lockdown modeling, you can implement detailed NPIs with location-specific dampings as in the SECIR model. The SECIRVVS model supports the same contact locations (e.g., home, school, work, other) and can apply different dampings to each location.

Example of defining locations and interventions for a detailed scenario:

.. code-block:: cpp

    // Define different contact locations
    enum class ContactLocation {
        Home = 0,
        School,
        Work,
        Other,
        Count,
    };
    
    // Define intervention types
    enum class Intervention {
        Home,
        SchoolClosure,
        HomeOffice,
        GatheringBanFacilitiesClosure,
        PhysicalDistanceAndMasks,
        SeniorAwareness,
    };

The model also supports dynamic NPIs based on epidemic thresholds:

.. code-block:: cpp

    // Configure dynamic NPIs
    auto& dynamic_npis = params.get<mio::osecirvvs::DynamicNPIsInfectedSymptoms<double>>();
    dynamic_npis.set_interval(mio::SimulationTime(3.0));  // Check every 3 days
    dynamic_npis.set_duration(mio::SimulationTime(14.0)); // Apply for 14 days
    dynamic_npis.set_base_value(100'000);                // Per 100,000 population
    dynamic_npis.set_threshold(200.0, dampings);         // Trigger at 200 cases per 100,000

Simulation
----------

The SECIRVVS model offers the same simulation functions as the SECIR model:

1. **simulate**: Standard simulation that tracks the compartment sizes over time
2. **simulate_flows**: Extended simulation that additionally tracks the flows between compartments

Basic simulation:

.. code-block:: cpp

    double t0 = 0;       // Start time
    double tmax = 30;    // End time
    double dt = 0.1;     // Time step
    
    // Run a standard simulation
    mio::TimeSeries<double> result = mio::osecirvvs::simulate<double>(t0, tmax, dt, model);

During simulation, the model handles several special processes:

1. **Daily Vaccinations**: At the beginning of each simulated day, the `apply_vaccination` function updates the immunity level of susceptible individuals based on the vaccination parameters.

2. **Variant Evolution**: If configured, the `apply_variant` function updates the transmission probability based on the existance of a new variant over time.

For both simulation types, you can also specify a custom integrator:

.. code-block:: cpp

    auto integrator = std::make_unique<mio::RKIntegratorCore>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    
    mio::TimeSeries<double> result = mio::osecirvvs::simulate(t0, tmax, dt, model, std::move(integrator));

Output
------

The output of the simulation is a `TimeSeries` object containing the sizes of each compartment at each time point. For a basic simulation, you can access the results as follows:

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
    
    // Print with custom column labels (you'll need a lot for this model!)
    std::vector<std::string> labels = {"S_naive", ... };
    result.print_table(std::cout, labels);

The order of the compartments is as defined in the `InfectionState` enum.

Additionally, you can export the results to a CSV file for further analysis or visualization:

.. code-block:: cpp

    // Export results to CSV
    result.export_csv("simulation_results.csv", labels);

Visualization
-------------

To visualize the results of a simulation, you can use the Python package :doc:`memilio_plot <../../python/memilio_plot>`
and its documentation. You can export your simulation results to CSV format as described above.

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
