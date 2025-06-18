Temporal-hybrid model
======================

Introduction
------------

The temporal-hybrid model switches between two models during the course of the simulation according to a given condition. Both models have to be initialized with their corresponding parameters and are handed to the class ``TemporalHybridSimulation``. 

The ``TemporalHybridSimulation`` class needs the used model types as well as their result types as template arguments. The results of the models are used to evaluate the switching condition. Additionally, conversion functions to convert the first model to the second model and vice versa have to be implemented for the used model types.
We implemented conversion functions for the following model combinations:

- Diffusive agent-based model, see ``mio::dabm``, using the singlewell potential and the stochastic metapopulation model, see ``mio::smm``.
- Diffusive agent-based model, see ``mio::dabm``, using the singlewell potential, and the ODE-based SECIR-type model, see ``mio::osecir``.

Simulation
----------

The simulation runs in discrete time steps. In every step, it is first checked whether the model needs to be switched by evaluating the switching condition. If the condition is fulfilled, the current model is converted into the target model i.e. its current state/result is used to update the state of the target model. Then the (new) currently used model is advanced until the next time step.

For a detailed description and application of the model, see:

- Bicker J, Schmieding R, et al. (2025) Hybrid metapopulation agent-based epidemiological models for efficient insight on the individual scale: A contribution to green computing. Infectious Disease Modelling, Volume 10, Issue 2. https://doi.org/10.1016/j.idm.2024.12.015

How to: Set up and run a simulation of the temporal-hybrid model
----------------------------------------------------------------

The following example shows how to run a temporal-hybrid simulation of the diffusive agent-based model using the singlewell potential and the ode-secir model.

.. code-block:: cpp

    using ABM = mio::dabm::Model<SingleWell<mio::hybrid::InfectionState>>;
    using ODE = mio::osecir::Model<double>;

First, both models have to be initialized. Please see the documentation of the corresponding models for information about their parameters and how to initialize them. In the following example code block, we initialize a ode-secir model with one age group and a diffusive ABM. We assume that the parameters of the diffusive ABM were created and initialized before and the parameters of the ode-secir model are set such that they match the ABM parameters.

.. code-block:: cpp

    //Initialize ABM
    ABM abm(agents, adoption_rates, interaction_radius, noise,
            {mio::hybrid::InfectionState::InfectedSevere, mio::hybrid::InfectionState::InfectedCritical,
             mio::hybrid::InfectionState::Dead});
    //Initialize ODE Model
    ODE ode(1);

After initializing the models, the corresponding simulations have to be created and the functions returning their result/current status have to be defined. For both models, the result function used for the hybridization just returns their result time series.

.. code-block:: cpp

    //Set t0 and internal dt for each model
    double t0 = 0;
    double dt = 0.1;

    //Create simulations
    auto sim_abm = mio::dabm::Simulation(abm, t0, dt);
    auto sim_ode = mio::Simulation(ode, t0, dt);

    const auto result_fct_abm = [](const mio::dabm::Simulation<SingleWell<mio::hybrid::InfectionState>>& sim,
                                   double /*t*/) {
        return sim.get_result();
    };

    const auto result_fct_ode = [](const mio::Simulation<double, ODE>& sim, double /*t*/) {
        return sim.get_result();
    };

Initializing the temporal-hybrid simulation with a given step size `dt_switch` with which the switching condition is checked is done as follows:

.. code-block:: cpp

    //Create hybrid simulation
    double dt_switch = 0.2;
    mio::hybrid::TemporalHybridSimulation<decltype(sim_abm), decltype(sim_ode), mio::TimeSeries<double>,
                                          mio::TimeSeries<double>>
        hybrid_sim(sim_abm, sim_ode, result_fct_abm, result_fct_ode, true, t0, dt_switch);

Before advancing the simulation until `tmax`, a switching condition has to be defined. In the example below, the temporal-hybrid model should switch from ABM to ODE if the number of infected individuals is bigger than 20 and it should switch back if the number is below 20.

.. code-block:: cpp

        //Define switching condition
    const auto condition = [](const mio::TimeSeries<double>& result_abm, const mio::TimeSeries<double>& result_ode,
                              bool abm_used) {
        if (abm_used) {
            auto& last_value = result_abm.get_last_value().eval();
            if ((last_value[(int)mio::hybrid::InfectionState::Exposed] +
                 last_value[(int)mio::hybrid::InfectionState::InfectedNoSymptoms] +
                 last_value[(int)mio::hybrid::InfectionState::InfectedSymptoms] +
                 last_value[(int)mio::hybrid::InfectionState::InfectedSevere] +
                 last_value[(int)mio::hybrid::InfectionState::InfectedCritical]) > 20) {
                return true;
            }
        }
        else {
            auto& last_value = result_ode.get_last_value().eval();
            if ((last_value[(int)mio::osecir::InfectionState::Exposed] +
                 last_value[(int)mio::osecir::InfectionState::InfectedNoSymptoms] +
                 last_value[(int)mio::osecir::InfectionState::InfectedNoSymptomsConfirmed] +
                 last_value[(int)mio::osecir::InfectionState::InfectedSymptoms] +
                 last_value[(int)mio::osecir::InfectionState::InfectedSymptomsConfirmed] +
                 last_value[(int)mio::osecir::InfectionState::InfectedSevere] +
                 last_value[(int)mio::osecir::InfectionState::InfectedCritical]) <= 20) {
                return true;
            }
        }
        return false;
    };

    //Simulate for 30 days
    double tmax = 30.;
    hybrid_sim.advance(tmax, condition);

The result ``mio::TimeSeries`` objects of the two models used (which are returned by the above defined result functions) can be accessed and printed via

.. code-block:: cpp

    //Print result time series of both models
    auto ts_abm = hybrid_sim.get_result_model1();
    auto ts_ode = hybrid_sim.get_result_model2();

    ts_abm.print_table({"S", "E", "Ins", "Isy", "Isev", "Icri", "R", "D"});
    ts_ode.print_table({"S", "E", "Ins", "Ins_confirmed", "Isy", "Isy_confirmed", "Isev", "Icri", "R", "D"});
