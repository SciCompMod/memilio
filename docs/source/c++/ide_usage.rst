IDE Models
==========

Models based on integro-differential equations allow arbitrary stay time distributions...

In MEmilio, two models based on integro-differential equations are implemented. They have different infection states and are solved with different numerical solvers. Their respective usage is described below.


IDE-SECIR model
----------------

This IDE-based model implements eight infection states and corresponding transition distributions can be set in a flexible way. The model is solved with a nonstandard finite difference scheme

For a detailed description and application of the model, see:

Wendler A, Plötzke L, Tritzschak H, Kühn MJ (2024) A nonstandard numerical scheme for a novel SECIR integro-differential equation-based model allowing nonexponentially distributed stay times. Under review. https://doi.org/10.48550/arXiv.2408.12228


How to: Set up and run a simulation of the IDE-SECIR model
-----------------------------------------------------------

The numerical solver requires a fixed time step size which we define by

.. code-block:: cpp

    ScalarType dt = 1.;

To initialize the model, the following inputs need to be passed to the model constructor:

- a time series containing the flows within a time step between the infection states for a large enough number of time points before the start of the simulation,
- a vector containing the population sizes for every age group,
- a vector containing the total number of deaths at time t0 for every age group,
- the number of age groups,
- optionally, a vector containing the total confirmed cases at time t0 and can be set if it should be used for initialization for every age group.

The number of age groups, the population sizes and total number of deaths can be defined directly by 

.. code-block:: cpp

    size_t num_agegroups = 1;

    mio::CustomIndexArray<ScalarType, mio::AgeGroup> N =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 10000.);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), 13.10462213);

In this example, we define the necessary flows before the simulation start by defining a time series at time points -10, ..., 0 that all contain the same vector of flows. The number of required time points before the simulation start depends on the chosen transition distributions that we can adapt later. 
Note that the last time point in our initial flow TimeSeries determines the start time of the simulation. 

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

With this, we can construct our model:

.. code-block:: cpp
    mio::isecir::Model model(std::move(init), N, deaths, num_agegroups);


TODO: Mention other init possibilities. 


If we do not want to use the default parameters, we can adapt them as follows. 

An important feature of our IDE-based model is that we can choose the transition distributions in a flexible way. The default distribution is a smoother cosine function as it provides good testing qualities. For more realistic simulations, MEmilio provides the possibility to use exponential, gamma or lognormal distributions within the model.
Practically, one first needs to create an object of a class that is derived from the class ``StateAgeFunction``, e.g. ``SmootherCosine``. Any class that is derived from ``StateAgeFunction`` can be inserted into a ``StateAgeFunctionWrapper`` object that is then passed to the model.

In this example, we start with creating a ``SmootherCosine`` object with parameter 2.0 that is then passed to the ``StateAgeFunctionWrapper`` object. Then we create a vector of type ``StateAgeFunctionWrapper``. Within this vector we adapt the distribution parameter for the transition from ``InfectedNoSymptoms`` to ``InfectedSymptoms``. Finally, this vector of ``StateAgeFunctionWrapper`` objects is passed to the model as demosntarted below.

.. code-block:: cpp
    mio::SmootherCosine smoothcos(2.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
        .set_distribution_parameter(4.0);

    model.parameters.get<mio::isecir::TransitionDistributions>()[mio::AgeGroup(0)] = vec_delaydistrib;

The transition probabilities can be set as follows

.. code-block:: cpp
    std::vector<ScalarType> vec_prob(num_transitions, 0.5);
    // The following probabilities must be 1, as there is no other way to go.
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
    model.parameters.get<mio::isecir::TransitionProbabilities>()[mio::AgeGroup(0)]        = vec_prob;

Here, we set the contact matrix used in the simulation. One can define multiple matrices for different locations. The size of each of these matrices is defined by tha number of age groups. 
In our example below we use only one contact matrix. We only consider one age group and set the contact rate to 10. 

.. code-block:: cpp
    size_t num_matrices =1;
    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(num_matrices, num_agegroups);
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_agegroups, num_agegroups, 10.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);


The parameters ``TransmissionProbabilityOnContact``, ``RelativeTransmissionNoSymptoms`` and ``RiskOfInfectionFromSymptomatic`` can be made dependent on the time. This is why we use a ``StateAgeFunction`` that is passed to a ``StateAgeFunctionWrapper`` to set these parameters. Note that there is also a ``ConstantFunction`` availbale if we do not want to have any dependency on the time. 
Here we use an ``ExponentialSurvivalFunction`` to set the mentioned parameters. 

.. code-block:: cpp
    mio::ExponentialSurvivalFunction exponential(0.5);
    mio::StateAgeFunctionWrapper prob(exponential);

    model.parameters.get<mio::isecir::TransmissionProbabilityOnContact>()[mio::AgeGroup(0)] = prob;
    model.parameters.get<mio::isecir::RelativeTransmissionNoSymptoms>()[mio::AgeGroup(0)]   = prob;
    model.parameters.get<mio::isecir::RiskOfInfectionFromSymptomatic>()[mio::AgeGroup(0)]   = prob;

Finally, we can also set the parameters Seasonality and StartDay directly as follows. 

.. code-block:: cpp
    model.parameters.set<mio::isecir::Seasonality>(0.1);
    // Start the simulation on the 40th day of a year (i.e. in February).
    model.parameters.set<mio::isecir::StartDay>(40);

Before the simulation, we check if all constraints of the model are satisfied so that the simulation can run as expected. 

.. code-block:: cpp
    model.check_constraints(dt);

To simulate the model from `t0` (that is determined by the initial flows provided to the constructor) to `tmax` with given step size `dt`, a Simulation has to be created and advanced until `tmax`, which is done as follows: 

.. code-block:: cpp
    ScalarType tmax = 10.;

    mio::isecir::Simulation sim(model, dt);
    sim.advance(tmax);

We can access and print the computed compartments and flows as follows. 

.. code-block:: cpp
    auto compartments = sim.get_result();
    auto flows = sim.get_transitions();

    compartments.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);
    flows.print_table({"S->E 1", "E->C 1", "C->I 1", "C->R 1", "I->H 1", "I->R 1", "H->U 1", "H->R 1", "U->D 1", "U->R 1"}, 16, 8);

If one wants to interpolate the results to a TimeSeries containing only full days, this can be done by

.. code-block:: cpp
    auto interpolated_results = mio::interpolate_simulation_result(sim.get_result());


TODO:   
- different initialization possible




IDE-SEIR model
---------------
This IDE-based model implements four infection states. 

For a detailed description and application of the model, see:

Ploetzke ... BA

Simulation
-----------

How to: Set up and run a simulation of the IDE-SEIR model
----------------------------------------------------------