LCT-based SECIR-type model
==========================

The LCT-SECIR module models and simulates an epidemic using an ODE-based approach while making use of the Linear Chain 
Trick to provide the option of Erlang distributed stay times in the compartments through the use of subcompartments. 
The model is particularly suited for pathogens with pre- or asymptomatic infection states and when severe or critical 
states are possible. The model assumes perfect immunity after recovery and is thus only suited for epidemic use cases. 
In the following, we present the model in detail.

Below is a visualization of the infection states and transitions. without a stratification according to groups.

.. image:: https://github.com/SciCompMod/memilio/assets/70579874/6a5d5a95-20f9-4176-8894-c091bd48bfb7
   :alt: tikzLCTSECIR


For a detailed description and application of the model, see:

- Plötzke L, Wendler A, Schmieding R, Kühn MJ. (2024). *Revisiting the Linear Chain Trick in epidemiological models: 
Implications of underlying assumptions for numerical solutions*. Under review. `https://doi.org/10.48550/arXiv.2412.09140 <https://doi.org/10.48550/arXiv.2412.09140>`_ 
- Plötzke L. (2023). *Der Linear Chain Trick in der epidemiologischen Modellierung als Kompromiss zwischen gewöhnlichen 
und Integro-Differentialgleichungen*. Master's thesis, University of Cologne. `https://elib.dlr.de/200381/ <https://elib.dlr.de/200381/>`_
- Hurtado PJ, Kirosingh AS. (2019). *Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell 
time distributions into mean field ODE models*. Journal of Mathematical Biology. `https://doi.org/10.1007/s00285-019-01412-w <https://doi.org/10.1007/s00285-019-01412-w>`_


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

It is possible to include subcompartments for the five compartments `Exposed`, `InfectedNoSymptoms`, `InfectedSymptoms`, `InfectedSevere` and `InfectedCritical`.


Sociodemographic Stratification
-------------------------------

In the LDE-SECIR model, the population can be stratified by one sociodemographic dimension. This dimension is denoted 
**Group**. It can be used for age groups as well as for other interpretations. 


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


Initial conditions
------------------

To initialize the model, we start by defining the number of subcompartments and constructing the model. We can choose the number of subcompartments.

.. code-block:: cpp
    
    constexpr size_t NumExposed = 2, NumInfectedNoSymptoms = 3, NumInfectedSymptoms = 1, NumInfectedSevere = 1,
                     NumInfectedCritical = 5;
    using InfState                       = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms,
                                            NumInfectedSevere, NumInfectedCritical, 1, 1>;
    using Model    = mio::lsecir::Model<LctState>;
    Model model;

For the simulation, we need initial values for all (sub)compartments. If we do not set the initial values manually, these are internally set to :math:`0`.

We start with constructing a vector ``initial_populations`` that we will pass on to the model. It contains vectors for each compartment, that contains a vector with initial values for the respective subcompartments. 
    
.. code-block:: cpp

        std::vector<std::vector<ScalarType>> initial_populations = {{750}, {30, 20},          {20, 10, 10}, {50},
                                                                    {50},  {10, 10, 5, 3, 2}, {20},         {10}};

We assert that vector has the correct size by checking that the number of `InfectionState`\s and the numbers of subcompartments are correct.

.. code-block:: cpp

        if (initial_populations.size() != (size_t)InfState::Count) {
            mio::log_error(
                "The number of vectors in initial_populations does not match the number of InfectionStates.");
            return 1;
        }
        if ((initial_populations[(size_t)InfState::Susceptible].size() !=
             LctState::get_num_subcompartments<InfState::Susceptible>()) ||
            (initial_populations[(size_t)InfState::Exposed].size() != NumExposed) ||
            (initial_populations[(size_t)InfState::InfectedNoSymptoms].size() != NumInfectedNoSymptoms) ||
            (initial_populations[(size_t)InfState::InfectedSymptoms].size() != NumInfectedSymptoms) ||
            (initial_populations[(size_t)InfState::InfectedSevere].size() != NumInfectedSevere) ||
            (initial_populations[(size_t)InfState::InfectedCritical].size() != NumInfectedCritical) ||
            (initial_populations[(size_t)InfState::Recovered].size() !=
             LctState::get_num_subcompartments<InfState::Recovered>()) ||
            (initial_populations[(size_t)InfState::Dead].size() !=
             LctState::get_num_subcompartments<InfState::Dead>())) {
            mio::log_error(
                "The length of at least one vector in initial_populations does not match the related number of "
                "subcompartments.");
            return 1;
        }

Now, we transfer the vector ``initial_populations`` to the model. 

.. code-block:: cpp

        std::vector<ScalarType> flat_initial_populations;
        for (auto&& vec : initial_populations) {
            flat_initial_populations.insert(flat_initial_populations.end(), vec.begin(), vec.end());
        }
        for (size_t i = 0; i < LctState::Count; i++) {
            model.populations[i] = flat_initial_populations[i];
        }
    }

    
In addition to setting the initial populations manually, MEmilio provides two other ways of setting the initial populations:  

- The file `parameters_io <https://github.com/SciCompMod/memilio/blob/main/cpp/models/lct_secir/parameters_io.h>`_ provides 
functionality to compute an initial value vector for the LCT-SECIR model based on reported data.
- The file `initializer_flows <https://github.com/SciCompMod/memilio/blob/main/cpp/models/lct_secir/initializer_flows.h>`_ 
provides functionality to compute an initial value vector for the LCT-SECIR model based on initial data in the form of 
a TimeSeries of InfectionTransitions. For the concept of the InfectionTransitions or flows, see also the IDE-SECIR model. 
This method can be particularly useful if a comparison is to be made with an IDE model with matching initialization or 
if the reported data is in the form of flows.


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


For more complex scenarios, such as real-world lockdown modeling, you can implement detailed NPIs with location-specific 
dampings. The SECIR model supports contact matrices for different locations (e.g., home, school, work, other) and can apply different dampings to each location.

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

We can simulate using the defined model from :math:`t_0` to :math:`t_{\max}` with initial step size :math:`dt` as follows:

.. code-block:: cpp

    ScalarType t0 = 0;
    ScalarType tmax = 10;
    ScalarType dt = 0.5;
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(t0, tmax, dt, model);

You can also specify a custom integrator:

.. code-block:: cpp

    auto integrator = std::make_shared<mio::RKIntegratorCore>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(t0, tmax, dt, model, integrator);


Output
------

The simulation result is divided by subcompartments. We can call the function ``calculate_compartments()`` to get a 
result according to the `InfectionState`\s .

.. code-block:: cpp

    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_compartments(result);

You can access the data in the `TimeSeries` object as follows:

.. code-block:: cpp

    // Get the number of time points.
    auto num_points = static_cast<size_t>(result.get_num_time_points());
    
    // Access data at a specific time point.
    Eigen::VectorX value_at_time_i = result.get_value(i);
    ScalarType time_i = result.get_time(i);
    
    // Access the last time point.
    Eigen::VectorX last_value = result.get_last_value();
    ScalarType last_time = result.get_last_time();


You can print the simulation results as a formatted table:

.. code-block:: cpp

    // Print results to console with default formatting.
    result.print_table();
    
    // Print with custom column labels.
    std::vector<std::string> labels = {"S", "E", "C", "I", "H", "U", "R", "D"};
    result.print_table(labels);

Additionally, you can export the results to a CSV file:

.. code-block:: cpp

    // Export results to CSV with default settings.
    result.export_csv("simulation_results.csv");


Visualization
-------------

To visualize the results of a simulation, you can use the Python package :doc:`memilio_plot <../../python/memilio_plot>`
and its documentation.

You can export your simulation results to CSV format as described above.

    
Examples
--------

An example can be found at:

- `examples/lct_secir.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/lct_secir.cpp>`_ 


Overview of the ``lsecir`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::lsecir