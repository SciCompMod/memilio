ODE-based SECIR-type model
===========================

The ODE-SECIR module models and simulates an epidemic using an ODE-based SECIR-type model approach.
The model is particularly suited for pathogens with pre- or asymptomatic infection states and when severe or critical
symptoms are possible. The model assumes perfect immunity after recovery. It is thus only suited for epidemic use cases
and, mostly, early epidemic phases. 

*   A generalization of the model that allows for Gamma or Erlang distributed stay times is the :doc:`LCT-SECIR model <cpp/lsecir>`.
*   A generalization of the model that allows for arbitrary distributed stay times is the :doc:`IDE-SECIR model <cpp/isecir>`.
*   A generalization of the model that includes three immunity layers and vaccination is the :doc:`ODE-SECIRVVS model <cpp/osecirvvs>`.
*   A generalization of the model that includes three immunity layers, vaccination, and waning immunity is the :doc:`ODE-SECIRTS model <cpp/osecirts>`.

The infection states and the transitions (also see next two sections) are visualized in the following graph.

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

The ODE-SECIR model is implemented as a **FlowModel**, which defines the derivatives of each flow between compartments.
This allows for explicit computation of new transmissions, infections, and hospitalizations. Additionally, the aggregated
compartment values can be computed with minimal overhead. The defined transitions `FromState, ToState` are:

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

The model implements the following parameters:

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
     - ``RiskOfInfectionFromSymptomatic``
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

The initial conditions of the model are represented by the class **Populations** which defines the number of individuals in each sociodemographic group and **InfectionState**. Before running a simulation, you need to set the initial values for each compartment:

.. code-block:: cpp

    // Set total population size
    model.populations.set_total(nb_total_t0); 
    
    // Set values for each InfectionState in the specific age group
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}] = nb_exp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] = nb_car_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] = nb_inf_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}] = nb_hosp_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}] = nb_icu_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}] = nb_rec_t0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}] = nb_dead_t0;
    
    // Set the susceptible population as difference to ensure correct total population
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}, nb_total_t0);

For age-resolved models, you need to set the initial conditions for each age group. Additionally, you can use ``set_difference_from_group_total`` to set the susceptible compartment as the difference between the total group size and all other compartments:

.. code-block:: cpp

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        model.populations[{i, mio::osecir::InfectionState::Exposed}] = fact * nb_exp_t0;
        // ...other states...
        model.populations.set_difference_from_group_total<mio::AgeGroup>(
            {i, mio::osecir::InfectionState::Susceptible}, fact * nb_total_t0);
    }


.. _Nonpharmaceutical Interventions:
Nonpharmaceutical Interventions
-------------------------------

In the SECIR model, nonpharmaceutical interventions (NPIs) are implemented through dampings to the contact matrix. These dampings reduce the contact rates between different groups to simulate interventions.

Basic dampings can be added to the contact matrix as follows:

.. code-block:: cpp

    // Create a contact matrix with constant contact rates between all groups
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    
    // Add a damping that reduces contacts by 70% starting at day 30
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

For age-resolved models, you can apply different dampings to different groups:

.. code-block:: cpp

    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    
    // Add a damping that reduces contacts within the same age group by 70% starting at day 30
    contact_matrix.add_damping(Eigen::VectorX<ScalarType>::Constant((size_t)nb_groups, 0.7).asDiagonal(),
                             mio::SimulationTime(30.));

The SECIR model also supports dynamic NPIs based on epidemic thresholds. These are implemented in the model specific **Simulation** class and are automatically triggered based on predefined criteria, such as the percentage of infected individuals in the population.

For more complex scenarios, such as real-world lockdown modeling, you can implement detailed NPIs with location-specific dampings. The SECIR model supports contact matrices for different locations (e.g., home, school, work, other) and can apply different dampings to each location.

Example for defining different contact locations:

.. code-block:: cpp

    // Define different contact locations
    enum class ContactLocation
    {
        Home = 0,
        School,
        Work,
        Other,
        Count,
    };
    
    // Map contact locations to strings for loading data files
    const std::map<ContactLocation, std::string> contact_locations = {
        {ContactLocation::Home, "home"},
        {ContactLocation::School, "school_pf_eig"},
        {ContactLocation::Work, "work"},
        {ContactLocation::Other, "other"}
    };

You can create intervention types that target specific locations with different intensities:

.. code-block:: cpp

    // Different types of NPI
    enum class Intervention
    {
        Home,
        SchoolClosure,
        HomeOffice,
        GatheringBanFacilitiesClosure,
        PhysicalDistanceAndMasks,
        SeniorAwareness,
    };
    
    // Different levels of NPI
    enum class InterventionLevel
    {
        Main,
        PhysicalDistanceAndMasks,
        SeniorAwareness,
        Holidays,
    };

A complex lockdown scenario with multiple interventions starting on a specific date can be implemented via:

.. code-block:: cpp

    auto start_lockdown_date = mio::Date(2020, 3, 18);
    auto start_lockdown = mio::SimulationTime(mio::get_offset_in_days(start_lockdown_date, start_date));
    
    // Apply different dampings for each intervention type
    contact_dampings.push_back(contacts_at_home(start_lockdown, 0.6, 0.8));
    contact_dampings.push_back(school_closure(start_lockdown, 1.0, 1.0));
    contact_dampings.push_back(home_office(start_lockdown, 0.2, 0.3));
    contact_dampings.push_back(social_events(start_lockdown, 0.6, 0.8));
    contact_dampings.push_back(physical_distancing(start_lockdown, 0.4, 0.6));

A more advances structure to automatically activate interventions based on threshold criteria is given by **DynamicNPIs**.
Dynamic NPIs can be configured to trigger when the number of symptomatic infected individuals exceeds a certain relative threshold in the population. 
In contrast to static NPIs which are active as long as no other NPI gets implemented, dynamic NPIs are checked at regular intervals and get 
activated for a defined duration when the threshold is exceeded. As above, different dampings `contact_dampings` can be assigned to different contact locations
and are then triggered all at once the threshold is exceeded.
The following example shows how to set up dynamic NPIs based on the number of 200 symptomatic infected individuals per 100,000 population. 
It will be active for at least 14 days and checked every 3 days. If the last check after day 14 is negative, the NPI will be deactivated.

.. code-block:: cpp

    // Configure dynamic NPIs with thresholds
    auto& dynamic_npis = params.get<mio::osecir::DynamicNPIsInfectedSymptoms<double>>();
    dynamic_npis.set_interval(mio::SimulationTime(3.0));  // Check every 3 days
    dynamic_npis.set_duration(mio::SimulationTime(14.0)); // Apply for 14 days
    dynamic_npis.set_base_value(100'000);                // Per 100,000 population
    dynamic_npis.set_threshold(200.0, contact_dampings);         // Trigger at 200 cases per 100,000


Simulation
----------

The SECIR model offers two simulation functions:

1. **simulate**: Standard simulation that tracks the compartment sizes over time
2. **simulate_flows**: Extended simulation that additionally tracks the flows between compartments

Standard simulation:

.. code-block:: cpp

    double t0 = 0;       // Start time
    double tmax = 50;    // End time
    double dt = 0.1;     // Time step
    
    // Run a standard simulation
    mio::TimeSeries<double> secir = mio::osecir::simulate(t0, tmax, dt, model);

Flow simulation for tracking transitions between compartments:

.. code-block:: cpp

    // Run a flow simulation to additionally track transitions between compartments
    auto result = mio::osecir::simulate_flows(t0, tmax, dt, model);
    // result[0] contains compartment sizes, result[1] contains flows

For both simulation types, you can also specify a custom integrator:

.. code-block:: cpp

    auto integrator = std::make_unique<mio::RKIntegratorCore>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    
    mio::TimeSeries<double> secir = mio::osecir::simulate(t0, tmax, dt, model, std::move(integrator));


Output
------

The output of the simulation is a `TimeSeries` object containing the sizes of each compartment at each time point. For a basic simulation, you can access the results as follows:

.. code-block:: cpp

    // Get the number of time points
    auto num_points = static_cast<size_t>(secir.get_num_time_points());
    
    // Access data at a specific time point
    Eigen::VectorXd value_at_time_i = secir.get_value(i);
    double time_i = secir.get_time(i);
    
    // Access the last time point
    Eigen::VectorXd last_value = secir.get_last_value();
    double last_time = secir.get_last_time();

For flow simulations, the result consists of two `mio::TimeSeries` objects, one for compartment sizes and one for flows:

.. code-block:: cpp

    auto result = mio::osecir::simulate_flows(t0, tmax, dt, model);
    
    // Access compartment sizes
    auto compartments = result[0];
    
    // Access flows between compartments
    auto flows = result[1];

You can print the simulation results as a formatted table:

.. code-block:: cpp

    // Print results to console with default formatting
    secir.print_table();
    
    // Print with custom column labels
    std::vector<std::string> labels = {"S", "E", "C", "C_confirmed", "I", "I_confirmed", "H", "U", "R", "D"};
    secir.print_table(labels);

Additionally, you can export the results to a CSV file:

.. code-block:: cpp

    // Export results to CSV with default settings
    secir.export_csv("simulation_results.csv");

The ODE-SECIR model also provides utility functions to extract specific measures, such as the reproduction number:

.. code-block:: cpp

    // Calculate R value at a specific time index
    auto r_at_index = mio::osecir::get_reproduction_number(time_idx, sim);
    
    // Calculate R values for the entire simulation
    Eigen::VectorXd r_values = mio::osecir::get_reproduction_numbers(sim);


Visualization
-------------

To visualize the results of a simulation, you can use the Python package :doc:`m-plot <../../python/m-plot>`
and its documentation.

    
Examples
--------

Different examples can be found at:

- `examples/ode_secir.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secir.cpp>`_
- `examples/ode_secir_ageres.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secir_ageres.cpp>`_
- `examples/ode_secir_parameter_study.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secir_parameter_study.cpp>`_  

Overview of the ``osecir`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::osecir
