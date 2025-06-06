IDE models
==========

MEmilio implements two models based on integro-differential equations (IDEs) with different infection states. IDE-based 
models are a generalization of ODE-based models. Whereas ODE-based models assume an exponential distribution regarding 
the time spent in an infection state, IDE-based models allow for arbitrary distributions.

In the following, we present the general structure of the IDE-SECIR model. 
The IDE-SEIR model follows a different structure, an introduction and a detailed example are provided in a separate 
section below. 

Infection states
----------------

Every model contains a list of **InfectionState**\s that define particular features of the subpopulations in the particular state.

.. code-block:: RST

    `State1`
    `State2`
    `...`


Infection state transitions
---------------------------

Additionally, we define **InfectionTransition**\s that define the possible transitions between the **InfectionState**\s.
When solving the model, we obtain results for the **InfectionTransition**\s as well.

.. code-block:: RST

    `Transition1`
    `Transition2`
    `...`


Sociodemographic stratification
-------------------------------

For the IDE-SECIR model, the population can also be stratified by one sociodemographic dimension. This dimension is denoted 
**AgeGroup** but can also be used for other interpretations. 

Parameters
----------

The parameters of the model are defined as structs and are combined in a class ``ParameterSet<Param1, Param2, ...>``.
We use different types of parameters to represent epidemiological parameters such as the distributions of the stay times in a 
compartment or the contact rates between different age groups. Most model parameters are constants that describe 
pathogen-specific characteristics (possibly resolved by sociodemographic groups) and are represented by a vector with a
value for each sociodemographic group. To model different contact rates between different sociodemographic groups, we
use a parameter denoted **ContactPatterns** of type **UncertainContactMatrix**. The **UncertainContactMatrix** contains
a set of contact matrices of arbitrary length and which can represent the different contact locations in the model like 
schools, workplaces, or homes. The matrices can be loaded or stored in the particular example.
In the **ContactPatterns**, each matrix element stores baseline contact rates :math:`c_{i,j}` between sociodemographic 
group :math:`i` to group :math:`j`. The dimension of the matrix is automatically defined by the model initiation and it is reduced 
to one value if no stratifcation is used. The values can be adjusted during the simulation, e.g., through implementing 
nonpharmaceutical interventions, see the section on :ref:`Nonpharmaceutical Interventions`.
Parameters can get accessed via ``model.parameters.get<Param<double>>()`` and set via either 
``model.parameters.get<Param<double>>() = value``. 


Initial conditions
------------------

The initial conditions consist of a **TimeSeries** that contains the flows per time step for some time interval before 
the simulation start. The length of this time interval depends on the chosen transition distributions and takes into 
account how long there is a relevant fraction of individuals remaining in the compartments. For more information, see 
the documentation of **StateAgeFunction**. Note that the last time point of the initial **TimeSeries** determines the 
start time of the simulation. 


.. _Nonpharmaceutical Interventions:
Nonpharmaceutical interventions
-------------------------------

Contact rates can be adjusted during the simulation to model nonpharmaceutical interventions (NPIs) such as lockdowns, 
school closures, or social distancing. This is done by adding **Damping**\s to the **ContactPatterns** of the model. A 
**Damping** is defined by a time point at which the intervention starts and a matrix of the same size as the 
**ContactMatrix**. While in many cases, the reduction matrix is given by a constant matrix with factor :math:`r`, also 
group-specific reductions are possible through setting particular rows or columns differently. With a constant 
reduction factor :math:`r`, the reduced contact rate is :math:`(1-r) * c_{i,j}`.

.. dropdown:: :fa:`gears` Expert's settings

    In some settings, contact rates cannot be reduced to zero to keep essential sectors of the society running. In this 
    case, we distinguish between a baseline contact matrix which we denote by :math:`B` and a minimum contact matrix 
    which we denote by :math:`M`. With a damping matrix :math:`D` the reduced contact matrix is then given by 
    :math:`B - D * (B - M)`, where the multiplication is element-wise.
    You can set the minimum and baseline contact matrices via 

    .. code-block:: cpp

        model.parameters.get<ContactPatterns>().get_cont_freq_mat()[0].get_baseline() = baseline_matrix;
        model.parameters.get<ContactPatterns>().get_cont_freq_mat()[0].get_minimum() = minimum_matrix;


Simulation
----------

Once the model is setup, run a simple simulation from time ``t0`` (determined by the initial conditions) to ``tmax`` 
with a fixed step size ``dt`` using the ``Simulation`` class defined in the ``models/ide_secir`` folder. The simulation 
uses a nonstandard numerical scheme to solve the IDEs that is implemented in MEmilio. 


Output
------

The output of the **Simulation** ``sim`` is a ``TimeSeries`` containing the sizes of each compartment at each time point 
and a ``TimeSeries`` containing the flows within a time step for each time point. A simple table can be printed using the 
``print_table()`` function of the ``TimeSeries`` class. The compartment sizes can be printed with 
``sim.get_result().print_table()`` and the flows with ``sim.get_transitions().print_table()``. 
As the time step may be small it may be useful to obtain outputs on days or user-defined time points. You can interpolate 
the results to days or any other series of times points with ``mio::interpolate_simulation_result()``.


Visualization
-------------

To visualize the results of a simulation, you can use the Python package :doc:`memilio_plot <../python/memilio_plot>`
and its documentation.    

List of models
-----------------------

.. toctree::
    :titlesonly:
    
    models/iseir
    models/isecir


In the following, we give detailed explanations both on the IDE-SECIR and the IDE-SEIR model. 

IDE-SECIR model
----------------

Introduction
~~~~~~~~~~~~~

This model is based on integro-differential equations.
The eight compartments

.. code-block::

    `Susceptible` (:math:`S`), may become Exposed at any time
    `Exposed` (:math:`E`), becomes InfectedNoSymptoms after some time
    `InfectedNoSymptoms` (:math:`I_{NS}`), becomes InfectedSymptoms or Recovered after some time
    `InfectedSymptoms` (:math:`I_{Sy}`), becomes InfectedSevere or Recovered after some time
    `InfectedSevere` (:math:`I_{Sev}`), becomes InfectedCritical or Recovered after some time
    `InfectedCritical` (:math:`I_{Cr}`), becomes Recovered or Dead after some time
    `Recovered` (:math:`R`)
    `Dead` (:math:`D`)

are used to simulate the spread of the disease and corresponding transition distributions can be set in a flexible way. 

Simulation
~~~~~~~~~~~

The simulation runs in discrete time steps using a nonstandard numerical scheme. This approach is based on the paper `"A non-standard numerical scheme for an age-of infection epidemic model" by Messina et al., Journal of Computational Dynamics, 2022 <https://doi.org/10.3934/jcd.2021029>`_. 

For a detailed description and application of the model, see:

Wendler A, Plötzke L, Tritzschak H, Kühn MJ (2024) A nonstandard numerical scheme for a novel SECIR integro-differential equation-based model allowing nonexponentially distributed stay times. Under review. https://doi.org/10.48550/arXiv.2408.12228


How to: Set up and run a simulation of the IDE-SECIR model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the following, we will demonstrate how to run a simulation using the IDE-SECIR model.

We start by defining the step size of the numerical solver. Note that the time step size is fixed during the simulation. Here, it is given by

.. code-block:: cpp

    ScalarType dt = 1.;

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

If we do not want to use the default parameters, we can adapt them as follows. 

An important feature of our IDE-based model is that we can choose the transition distributions in a flexible way. The default distribution is a smoother cosine function as it provides good testing qualities. For more realistic simulations, MEmilio provides the possibility to use exponential, gamma or lognormal distributions within the model.
Practically, one first needs to create an object of a class that is derived from the class ``StateAgeFunction``, e.g. ``SmootherCosine``. Any class that is derived from ``StateAgeFunction`` can be inserted into a ``StateAgeFunctionWrapper`` object that is then passed to the model.

In this example, we start with creating a ``SmootherCosine`` object with a parameter of :math:`2` that is then passed to the ``StateAgeFunctionWrapper`` object. Then we create a vector of type ``StateAgeFunctionWrapper``. Within this vector we adapt the distribution parameter for the transition from ``InfectedNoSymptoms`` to ``InfectedSymptoms``. Finally, this vector of ``StateAgeFunctionWrapper`` objects is passed to the model as demonstrated below.

.. code-block:: cpp

    mio::SmootherCosine smoothcos(2.0);
    mio::StateAgeFunctionWrapper delaydistribution(smoothcos);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
        .set_distribution_parameter(4.0);

    model.parameters.get<mio::isecir::TransitionDistributions>()[mio::AgeGroup(0)] = vec_delaydistrib;

The transition probabilities can be set as follows.

.. code-block:: cpp

    std::vector<ScalarType> vec_prob(num_transitions, 0.5);
    // The following probabilities must be 1, as there is no other way to go.
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
    model.parameters.get<mio::isecir::TransitionProbabilities>()[mio::AgeGroup(0)]        = vec_prob;

Now, we set the contact matrix determining the contacts between the age groups. One can define multiple matrices for different locations. The size of each of these matrices is defined by the number of age groups. 
Below, we use only one contact matrix for one location. As we only consider one age group in our example, we set the corresponding contact rate to :math:`10`. 

.. code-block:: cpp

    size_t num_matrices = 1;
    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(num_matrices, num_agegroups);
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_agegroups, num_agegroups, 10.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

The parameters ``TransmissionProbabilityOnContact``, ``RelativeTransmissionNoSymptoms`` and ``RiskOfInfectionFromSymptomatic`` can be made dependent on the time. This is why we use a ``StateAgeFunction`` that is passed to a ``StateAgeFunctionWrapper`` to set these parameters. Note that there is also a ``ConstantFunction`` available if we do not want to have any dependency on the time. 
Here we use an ``ExponentialSurvivalFunction`` to set the mentioned parameters. 

.. code-block:: cpp

    mio::ExponentialSurvivalFunction exponential(0.5);
    mio::StateAgeFunctionWrapper prob(exponential);

    model.parameters.get<mio::isecir::TransmissionProbabilityOnContact>()[mio::AgeGroup(0)] = prob;
    model.parameters.get<mio::isecir::RelativeTransmissionNoSymptoms>()[mio::AgeGroup(0)]   = prob;
    model.parameters.get<mio::isecir::RiskOfInfectionFromSymptomatic>()[mio::AgeGroup(0)]   = prob;

Finally, we can also set the parameters ``Seasonality`` and ``StartDay`` directly as follows. 

.. code-block:: cpp

    model.parameters.set<mio::isecir::Seasonality>(0.1);
    // Start the simulation on the 40th day of a year (i.e. in February).
    model.parameters.set<mio::isecir::StartDay>(40);

Before the simulation, we check if all constraints of the model are satisfied so that the simulation can run as expected. 

.. code-block:: cpp

    model.check_constraints(dt);

To simulate the model from :math:`t_0` (that is determined by the initial flows provided to the constructor) to :math:`t_{\max}` with given step size :math:`dt`, a object of the ``Simulation`` class has to be created and advanced until :math:`t_{\max}`, which is done as follows.

.. code-block:: cpp

    ScalarType tmax = 10.;

    mio::isecir::Simulation sim(model, dt);
    sim.advance(tmax);

We can access and print the computed compartments and flows. 

.. code-block:: cpp

    auto compartments = sim.get_result();
    auto flows = sim.get_transitions();

    compartments.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);
    flows.print_table({"S->E 1", "E->C 1", "C->I 1", "C->R 1", "I->H 1", "I->R 1", "H->U 1", "H->R 1", "U->D 1", "U->R 1"}, 16, 8);

If one wants to interpolate the results to a ``TimeSeries`` containing only full days, this can be done by

.. code-block:: cpp

    auto interpolated_results = mio::interpolate_simulation_result(sim.get_result());


IDE-SEIR model
---------------

Introduction
~~~~~~~~~~~~~
The four compartments 

- `Susceptible` (:math:`S`), may become exposed at any time
- `Exposed` (:math:`E`), becomes infected after some time
- `Infected` (:math:`I`), will recover after some time
- `Recovered` (:math:`R`)

are used to simulate the spread of the disease. 

Simulation
~~~~~~~~~~~

The simulation runs in discrete time steps using a trapezoidal rule. The model and the numerical scheme is based on the paper `"Modeling infectious diseases using integro-differential equations: Optimal
control strategies for policy decisions and Applications in COVID-19" by Keimer and Pflug, 2020 <http://dx.doi.org/10.13140/RG.2.2.10845.44000>`_. 

For a detailed description and application of the model, see:

Plötzke L (2021) Modellierung epidemischer Infektionskrankheiten auf der Basis von gewöhnlichen und Integro-Differentialgleichungen. Bachelor thesis, University of Cologne. https://elib.dlr.de/143504/

How to: Set up and run a simulation of the IDE-SEIR model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To initialize the model, the following inputs need to be passed to the model constructor:

- a ``TimeSeries`` containing the number of `Susceptible`s for a large enough number of time points before the start of the simulation,
- the time step size :math:`dt` used for numerical integration,
- the size of the population of the considered region :math:`N`. 

The initialization of the model can be done as follows where we set the `Susceptible`s from :math:`-15, \dots, 0` based on the total population and the time of the previous time point.

.. code-block:: cpp

    using Vec = mio::TimeSeries<double>::Vector;


    int N     = 810000;
    double dt = 0.1;
    mio::TimeSeries<double> init(1);

    init.add_time_point<Eigen::VectorXd>(-15.0, Vec::Constant(1, N * 0.95));
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt,
                            Vec::Constant(1, (double)init.get_last_value()[0] + init.get_last_time()));
    }

    // Initialize model.
    mio::iseir::Model<double> model(std::move(init), dt, N);

If we do not want to use the default parameters, we can adapt them as follows.

The parameters ``LatencyTime``, ``InfectiousTime`` and ``TransmissionRisk`` can be set directly. 

.. code-block:: cpp

    model.parameters.set<mio::iseir::LatencyTime>(3.3);
    model.parameters.set<mio::iseir::InfectiousTime>(8.2);
    model.parameters.set<mio::iseir::TransmissionRisk>(0.015);

Here, we set the contact matrix used in the simulation. One can define multiple matrices for different locations. The size of each of these matrices is defined by the number of age groups. 
Below, we use only one contact matrix for one location. As we only consider one age group in our example, we set the corresponding contact rate to :math:`10`.

.. code-block:: cpp

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));

To simulate the implementation of nonpharmaceutical interventions, we add dampings to the contact rate. Here, we apply a damping of :math:`0.7` after :math:`10` days, meaning that the contact rate is reduced to :math:`30%` of the initial value.  

.. code-block:: cpp

    contact_matrix[0].add_damping(0.7, mio::SimulationTime(10.));
    model.parameters.get<mio::iseir::ContactFrequency<double>>() = mio::UncertainContactMatrix<double>(contact_matrix);

After defining :math:`t_{\max}`, we can simulate, which means that we calculate the value for the compartment :math:`S`.

.. code-block:: cpp

    int tmax  = 15;
    model.simulate(tmax);

The values of the remaining compartments :math:`E`, :math:`I` and :math:`R` are calculated using the parameters ``LatencyTime`` and ``InfectiousTime`` and obtain a time series containing the values of all compartments. 

.. code-block:: cpp

    auto result = model.calculate_EIR();

Finally, we can print our results. 

.. code-block:: cpp

    result.print_table({"S", "E", "I", "R"});


