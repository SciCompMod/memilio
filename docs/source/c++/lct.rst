Linear Chain Trick model
=========================

Introduction
-------------

The Linear Chain Trick provides the option to use Erlang-distributed stay times in the compartments through the use of subcompartments. 
The normal ODE models have (possibly unrealistic) exponentially distributed stay times.
The LCT model can still be described by an ordinary differential equation system.

For a detailed description and application of the model, see:

- Plötzke L, Wendler A, Schmieding R, Kühn MJ (2024) Revisiting the Linear Chain Trick in epidemiological models: Implications of underlying assumptions for numerical solutions. Under review. https://doi.org/10.48550/arXiv.2412.09140
- Hurtado PJ und Kirosingh AS (2019) Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models. Journal of Mathematical Biology. https://doi.org/10.1007/s00285-019-01412-w

The eight compartments 

- `Susceptible` (:math:`S`), may become Exposed at any time
- `Exposed` (:math:`E`), becomes InfectedNoSymptoms after some time
- `InfectedNoSymptoms` (:math:`I_{NS}`), becomes InfectedSymptoms or Recovered after some time
- `InfectedSymptoms` (:math:`I_{Sy}`), becomes InfectedSevere or Recovered after some time
- `InfectedSevere` (:math:`I_{Sev}`), becomes InfectedCritical or Recovered after some time
- `InfectedCritical` (:math:`I_{Cr}`), becomes Recovered or Dead after some time
- `Recovered` (:math:`R`)
- `Dead` (:math:`D`)

are used to simulate the spread of the disease. 
It is possible to include subcompartments for the five compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.
You can divide the population according to different groups, e.g. AgeGroups or gender and choose parameters according to groups.

Simulation
-----------

How to: Set up and run a simulation of the LCT-SECIR model
-----------------------------------------------------------

We start by defining the number of subcompartments and constructing the model. 

.. code-block:: cpp
    
    // Simple example to demonstrate how to run a simulation using an LCT-SECIR model.
    // One single AgeGroup/Category member is used here.
    // Parameters, initial values and the number of subcompartments are not meant to represent a realistic scenario.
    constexpr size_t NumExposed = 2, NumInfectedNoSymptoms = 3, NumInfectedSymptoms = 1, NumInfectedSevere = 1,
                     NumInfectedCritical = 5;
    using InfState                       = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms,
                                            NumInfectedSevere, NumInfectedCritical, 1, 1>;
    using Model    = mio::lsecir::Model<LctState>;
    Model model;

    // Variable defines whether the class Initializer is used to define an initial vector from flows or whether a manually
    // defined initial vector is used to initialize the LCT model.
    bool use_initializer_flows = false;

    ScalarType tmax = 10;

    // Set Parameters.
    model.parameters.get<mio::lsecir::TimeExposed>()[0]            = 3.2;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0] = 2.;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]   = 5.8;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]     = 9.5;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]   = 7.1;

    model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[0] = 0.05;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    // From SimulationTime 5, the contact pattern is reduced to 30% of the initial value.
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.));

    model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[0] = 0.7;
    model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[0] = 0.25;
    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = 0.09;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = 0.2;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = 0.25;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = 0.3;

    if (use_initializer_flows) {
        // Example how to use the class Initializer for the definition of an initial vector for the LCT model.

        ScalarType dt                                    = 0.001;
        Eigen::VectorX<ScalarType> total_population      = Eigen::VectorX<ScalarType>::Constant(1, 1000000.);
        Eigen::VectorX<ScalarType> deaths                = Eigen::VectorX<ScalarType>::Constant(1, 10.);
        Eigen::VectorX<ScalarType> total_confirmed_cases = Eigen::VectorX<ScalarType>::Constant(1, 16000.);

        // Create TimeSeries with num_transitions elements.
        int num_transitions = (int)mio::lsecir::InfectionTransition::Count;
        mio::TimeSeries<ScalarType> flows(num_transitions);

        mio::TimeSeries<ScalarType>::Vector vec_flows(num_transitions);
        vec_flows[(int)mio::lsecir::InfectionTransition::SusceptibleToExposed]                 = 2.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 1.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 8.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 4.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 1.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedSymptomsToRecovered]          = 4.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 1.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedSevereToRecovered]            = 1.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedCriticalToDead]               = 1.0;
        vec_flows[(int)mio::lsecir::InfectionTransition::InfectedCriticalToRecovered]          = 1.0;
        vec_flows                                                                              = vec_flows * dt;
        // Add initial time point to time series.
        flows.add_time_point(-110, vec_flows);
        // Add further time points until time 0.
        while (flows.get_last_time() < -dt / 2) {
            flows.add_time_point(flows.get_last_time() + dt, vec_flows);
        }

        // Set initialization vector for the LCT model.
        mio::lsecir::Initializer<Model> initializer(std::move(flows), model);
        initializer.set_tol_for_support_max(1e-6);
        auto status = initializer.compute_initialization_vector(total_population, deaths, total_confirmed_cases);
        if (status) {
            return 1;
        }
    }
    else {
        // Simple example how to initialize model without flows.
        // Define the initial values with the distribution of the population into subcompartments.
        // This method of defining the initial values using a vector of vectors is not necessary, but should remind you
        // how the entries of the initial value vector relate to the defined template parameters of the model or the number
        // of subcompartments. It is also possible to define the initial values directly.
        std::vector<std::vector<ScalarType>> initial_populations = {{750}, {30, 20},          {20, 10, 10}, {50},
                                                                    {50},  {10, 10, 5, 3, 2}, {20},         {10}};

        // Assert that initial_populations has the right shape.
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

        // Transfer the initial values in initial_populations to the model.
        std::vector<ScalarType> flat_initial_populations;
        for (auto&& vec : initial_populations) {
            flat_initial_populations.insert(flat_initial_populations.end(), vec.begin(), vec.end());
        }
        for (size_t i = 0; i < LctState::Count; i++) {
            model.populations[i] = flat_initial_populations[i];
        }
    }

    // Perform a simulation.
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(0, tmax, 0.5, model);
    // The simulation result is divided by subcompartments.
    // We call the function calculate_compartments to get a result according to the InfectionStates.
    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_compartments(result);
    auto interpolated_results = mio::interpolate_simulation_result(population_no_subcompartments);
    interpolated_results.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 12, 4);




