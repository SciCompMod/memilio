.. include:: ../../literature.rst

IDE-based SECIR-type model
==========================

The IDE-SECIR module models and simulates an epidemic using integro-differential equations allowing for 
arbitrary stay time distributions in the compartments. The model is particularly suited for pathogens with pre- or 
asymptomatic infection states and when severe or critical states are possible. The model assumes perfect immunity after 
recovery. It is thus only suited for epidemic use cases and, mostly, early epidemic phases. 

*   A generalization of the model in the application sense that includes three immunity layers and vaccination is the :doc:`ODE-SECIRVVS model <cpp/osecirvvs>`. However, this model does not allow for arbitrary stay time distributions.
*   A generalization of the model in the application sense that includes three immunity layers, vaccination, and waning immunity is the :doc:`ODE-SECIRTS model <cpp/osecirts>`. However, this model does not allow for arbitrary stay time distributions.

Below is a visualization of the infection states and transitions. The variables :math:`\sigma_{z_1}^{z_2}` refer to a transition from a compartment :math:`z_1` to a compartment :math:`z_2`.

.. image:: https://github.com/SciCompMod/memilio/assets/70579874/3500421a-035c-4ce1-ae95-a54d8097be82
   :alt: tikzIDESECIR


The simulation runs in discrete time steps using a non-standard numerical scheme. This approach is based on the paper
Messina E, Pezzella M. (2022). *A non-standard numerical scheme for an age-of-infection epidemic model*. Journal of Computational Dynamics.
`https://doi.org/10.3934/jcd.2021029 <https://doi.org/10.3934/jcd.2021029>`_

A detailed investigation of the IDE-SECIR model and numerical experiments can be found in |A_nonstandard_numerical_scheme|

Infection States
----------------

The model contains the following list of **InfectionState**\s:

.. code-block:: RST

    `Susceptible`
    `Exposed`
    `InfectedNoSymptoms`
    `InfectedSymptoms`
    `InfectedSevere`
    `InfectedCritical`
    `Recovered`
    `Dead`

Infection State Transitions
---------------------------

The possible transitions between the **InfectionState**\s are:

.. code-block:: RST
  
    `SusceptibleToExposed`
    `ExposedToInfectedNoSymptoms`
    `InfectedNoSymptomsToInfectedSymptoms` 
    `InfectedNoSymptomsToRecovered`
    `InfectedSymptomsToInfectedSevere`
    `InfectedSymptomsToRecovered`
    `InfectedSevereToInfectedCritical`
    `InfectedSevereToRecovered`
    `InfectedCriticalToDead`
    `InfectedCriticalToRecovered`


Sociodemographic Stratification
-------------------------------

In the IDE-SECIR model, the population can be stratified by one sociodemographic dimension. This dimension is denoted 
**AgeGroup** but can also be used for other interpretations. 


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
     - Transmission risk for people located the susceptible compartment.
   * - :math:`\xi_{I_{NS}}`
     - ``RelativeTransmissionNoSymptoms``
     - Proportion of nonsymptomatically infected people who are not isolated.
   * - :math:`\xi_{I_{Sy}}`
     - ``RiskOfInfectionFromSymptomatic``
     - Proportion of infected people with symptoms who are not isolated.
   * - :math:`k`
     - ``Seasonality``
     - Seasonal variation factor affecting transmission.
   * - :math:`t_{start}`
     - ``StartDay``
     - Day of the year on which the simulation is started.
   * - :math:`\mu_{z_1}^{z_2}`
     - ``TransitionProbabilities``
     - Probability of transitioning from compartment :math:`z_1` to compartment :math:`z_2`.
   * - :math:`\gamma_{z_1}^{z_2}(\tau)`
     - ``TransitionDistributions``
     - Expected proportion of people who are still in compartment :math:`z_1` :math:`\tau` days after entering this compartment and who will move to compartment :math:`z_2` later in the course of the disease.


Initial conditions
------------------

To initialize the model, the following inputs need to be passed to the model constructor:

- a ``TimeSeries`` containing the flows within a time step between the infection states for a large enough number of time points before the start of the simulation,
- a vector containing the population sizes for every age group,
- a vector containing the total number of deaths at time :math:`t_0` for every age group,
- the number of age groups,
- optionally, a vector containing the total confirmed cases at time :math:`t_0` and can be set if it should be used for initialization for every age group.

The number of age groups, the population sizes and total number of deaths can be defined directly by 

.. code-block:: cpp

    size_t num_agegroups = 1;

    mio::CustomIndexArray<ScalarType, mio::AgeGroup> N =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 10000.);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 13.10462213);

In this example, we define the necessary flows before the simulation start by defining a time series at time points :math:`-10,\dots, 0` that all contain the same vector of flows. The number of required time points before the simulation start depends on the chosen transition distributions that we can adapt later. 
Note that the last time point in our initial flow ``TimeSeries`` determines the start time of the simulation. 

Note that the initial flows already require knowledge of the time step size of the numerical solver. This is foxed during the the simulation and can be set by

.. code-block:: cpp

    ScalarType dt = 1.;

Then we can define the initial flows as follows. 

.. code-block:: cpp

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // Create TimeSeries with num_transitions * num_agegroups elements where transitions needed for simulation will be
    // stored.
    mio::TimeSeries<ScalarType> init(num_transitions * num_agegroups);

    // Define vector with flows. 
    Vec vec_init(num_transitions * num_agegroups);
    vec_init[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 25.0;
    vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 15.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 4.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 1.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 4.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 1.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered]            = 1.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead]               = 1.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered]          = 1.0;

    vec_init = vec_init * dt;

    // Add initial time point to time series.
    init.add_time_point(-10, vec_init);
    // Add further time points until time 0.
    while (init.get_last_time() < -dt / 2) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

There are different options for initializing a fictional scenario. Regardless of the approach, you must provide a history of values for the transitions as demonstrated above and possibly additional information to compute the initial distribution of the population in the compartments. This information must be of the following type:  

    - You can state the number of total confirmed cases `total_confirmed_cases` at time :math:`t_0`. The number of recovered people is set accordingly and the remaining values are derived in the model before starting the simulation. Then the model can be constructed by 

    .. code-block:: cpp

        mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_confirmed_cases =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 100.);
        mio::isecir::Model model(std::move(init), N, deaths, num_agegroups, total_confirmed_cases);
    
    - If you cannot provide this number of total confirmed cases, we can construct the model without this information.

    .. code-block:: cpp
    
        mio::isecir::Model model(std::move(init), N, deaths, num_agegroups);

    In that case, we have three possible options for initializing:

        - You can set the number of people in the `Susceptible` compartment at time :math:`t_0` via `populations`. Initial values of the other compartments are derived in the model before starting the simulation.

        .. code-block:: cpp

            model.populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Susceptible] = 1000.;

        - You can set the number of people in the `Recovered` compartment at time :math:`t_0` via `populations`. Initial values of the other compartments are derived in the model before starting the simulation.

        .. code-block:: cpp

            model.populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Recovered] = 1000.;

        - If none of the above is used, the force of infection formula and the values for the initial transitions are used consistently with the numerical scheme proposed in `Messina et al (2022) <https://doi.org/10.3934/jcd.2021029>`_ to set the `Susceptible`s. 

- The file `parameters_io <https://github.com/SciCompMod/memilio/blob/main/cpp/models/ide_secir/parameters_io.h>`_ provides functionality to compute initial data for the IDE-SECIR model based on real data. An example for this initialization method can be found at  `IDE initialization example <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ide_initialization.cpp>`_.


.. _Nonpharmaceutical Interventions:
Nonpharmaceutical Interventions
-------------------------------

In the SECIR model, nonpharmaceutical interventions (NPIs) are implemented through dampings in the contact matrix. 
These dampings reduce the contact rates between different groups to simulate interventions.

Basic dampings can be added to the contact matrix as follows:

.. code-block:: cpp

    // Create a contact matrix with constant contact rates between all groups.
    ScalarType cont_freq = 10.;
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    
    // Add a uniform damping across all age groups.
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

For age-resolved models, you can apply different dampings to different groups:

.. code-block:: cpp

    ScalarType cont_freq = 10.;
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_agegroups, num_agegroups, cont_freq));
    
    // Add a damping that reduces contacts within the same age group by 70% starting at day 30.
    contact_matrix.add_damping(Eigen::VectorX<ScalarType>::Constant(num_agegroups, 0.7).asDiagonal(),
                             mio::SimulationTime(30.));


For more complex scenarios, such as real-world venue closures or lockdown modeling, you can implement detailed NPIs with location-specific dampings. The IDE-SECIR model supports contact matrices for different locations (e.g., home, school, work, other) and can apply different dampings to each location.

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


Simulation
----------

Before the simulation, we check if all constraints of the model are satisfied so that the simulation can run as expected. 

.. code-block:: cpp

    model.check_constraints(dt);

To simulate the model from :math:`t_0` (that is determined by the initial flows provided to the constructor) to 
:math:`t_{\max}` with given step size :math:`dt`, a object of the **Simulation** class has to be created and advanced 
until :math:`t_{\max}`, which is done as follows.

.. code-block:: cpp

    ScalarType tmax = 10.;

    mio::isecir::Simulation sim(model, dt);
    sim.advance(tmax);


Output
------

The output of the simulation are two `TimeSeries` objects, one containing the size of the compartments at all time 
points and one containing the number of transitions within a time step. You can access the results as follows:

.. code-block:: cpp
    
    // Access compartment sizes.
    auto compartments = sim.get_result();
    
    // Access transitions between compartments.
    auto transitions = sim.get_transitions();

The order of the compartments and transitions follows the definition in the **InfectionState** and **InfectionTransition** enums, respectively.

You can access the data in the `TimeSeries` objects as follows:

.. code-block:: cpp

    // Get the number of time points.
    auto num_points = static_cast<size_t>(compartments.get_num_time_points());
    
    // Access data at a specific time point.
    Eigen::VectorX value_at_time_i = compartments.get_value(i);
    ScalarType time_i = compartments.get_time(i);
    
    // Access the last time point.
    Eigen::VectorX last_value = compartments.get_last_value();
    ScalarType last_time = compartments.get_last_time();


You can print the simulation results as a formatted table:

.. code-block:: cpp

    // Print results to console with default formatting.
    compartments.print_table();
    
    // Print with custom column labels.
    std::vector<std::string> labels = {"S", "E", "C", "I", "H", "U", "R", "D"};
    compartments.print_table(labels);

Additionally, you can export the results to a CSV file:

.. code-block:: cpp

    // Export results to CSV with default settings.
    compartments.export_csv("simulation_results.csv");


Visualization
-------------

To visualize the results of a simulation, you can use the Python package :doc:`m-plot <../../python/m-plot>`
and its documentation.

You can export your simulation results to CSV format as described above.

    
Examples
--------

Different examples can be found at:

- `examples/ide_secir.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ide_secir.cpp>`_
- `examples/ide_secir_ageres.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ide_secir_ageres.cpp>`_
- `examples/ide_initialization.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ide_initialization.cpp>`_  

Overview of the ``isecir`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::isecir
