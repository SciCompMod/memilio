Linear Chain Trick model
=========================

MEmilio implements a SECIR-type model utilizing the Linear Chain Trick (LCT). This is a generalization of simple ODE 
models and allows for Erlang distributed stay times in the compartments by introducing subcompartments. Note that the 
resulting system can still be described by ODEs. The LCT-SECIR model can be stratified by age groups or other sociodemographic 
factors.

In the following, we present the general structure of the LCT model. The particular model documentation with examples 
is linked at the bottom of this page.

Infection states
----------------

The model contains a list of **InfectionState**\s that define particular features of the subpopulations in the particular state.

.. code-block:: RST

    `State1`
    `State2`
    `...`

To make use of the LCT, we additionally need to define the numbers of subcompartments for each **InfectionState**.

.. code-block:: RST
    `Number of subcompartments of State1`
    `Number of subcompartments of State2`
    `...`

Note that the model is implemented as **CompartmentalModel**.



Sociodemographic stratification
-------------------------------

For this model, the population can also be stratified by one sociodemographic dimension. This dimension can be used 
for age groups but can also be used for other interpretations. 


Parameters
----------

The parameters of the model are defined as structs and are combined in a class ``ParameterSet<Param1, Param2, ...>``.
We use different types of parameters to represent epidemiological parameters such as the mean stay times in a 
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
``model.parameters.get<Param<double>>() = value`` or ``model.parameters.set<Param<double>>(value)``. 


Initial conditions
------------------

The initial conditions of the model are represented by a class **LctPopulations** that gives the number of individuals in 
each sociodemographic group and each subcompartment of each **InfectionState**. For more details, see 
:doc:`Model Creation <lct_creation>`. Before the simulation, set the initial conditions for each **InfectionState** and 
sociodemographic group.


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

Once the model is setup, run a simple simulation from time ``t0`` to ``tmax`` with an initial step size ``dt`` using the 
``mio::simulation()`` function. This will run a simulation of type **Simulation** that saves the sizes of each 
subcompartment over time. 
You can run a simulation using either fixed or adaptive solution schemes with an absolute or relative tolerance. By 
default, the simulation uses an adaptive solution scheme of the boost library and an absolute tolerance of 1e-10 and a 
relative tolerance of 1e-5. For more details on the possible integration schemes, see <numerical integrator stuff>.


Output
------

The output of the **Simulation** is a **TimeSeries** ``result`` containing the sizes of each subcompartment at each time point. 
To obtain a result with respect to the compartments, the subcompartments can be accumulated via the function 
``calculate_compartments()``. A simple table can be printed using the ``print_table()`` function of the 
``TimeSeries`` class. The compartment sizes can be printed with ``model.calculate_compartments(result).print_table()``. 
As adaptive step size methods are used by default, the output will not be available on equidistant time points like `dt`
or days. To obtain outputs on days or user-defined time points, you can interpolate the results to days or
any other series of times points with ``mio::interpolate_simulation_result()``.


Visualization
-------------

To visualize the results of a simulation, you can use the Python package :doc:`memilio_plot <../python/memilio_plot>`
and its documentation.


In the following, we give detailed explanations of the LCT-SECIR model.

Introduction
-------------

The Linear Chain Trick (LCT) provides the option to use Erlang distributed stay times in the compartments through the 
use of subcompartments whereas simple ODE models have (possibly unrealistic) exponentially distributed stay times.
Note that LCT-based models can still be described by an ordinary differential equation system.

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
It is possible to include subcompartments for the five compartments `Exposed`, `InfectedNoSymptoms`, `InfectedSymptoms`, `InfectedSevere` and `InfectedCritical`.
You can divide the population according to different groups, e.g. by age or gender and choose parameters according to groups.

Simulation
-----------

The simulation runs in discrete time steps. Different ODE solvers are available, some of them use an adaptive time step size.

For a detailed description and application of the model, see:

- Plötzke L, Wendler A, Schmieding R, Kühn MJ (2024) Revisiting the Linear Chain Trick in epidemiological models: Implications of underlying assumptions for numerical solutions. Under review. https://doi.org/10.48550/arXiv.2412.09140
- Hurtado PJ und Kirosingh AS (2019) Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models. Journal of Mathematical Biology. https://doi.org/10.1007/s00285-019-01412-w

How to: Set up and run a simulation of the LCT-SECIR model
-----------------------------------------------------------

In the following, we will demonstrate how to run a simulation using the LCT-SECIR model. This examples uses one age group/category.

We start by defining the number of subcompartments and constructing the model. We can choose the number of subcompartments individually for the compartments `Exposed`, `InfectedNoSymptoms`, `InfectedSymptoms`, `InfectedSevere` and `InfectedCritical`.

.. code-block:: cpp
    
    constexpr size_t NumExposed = 2, NumInfectedNoSymptoms = 3, NumInfectedSymptoms = 1, NumInfectedSevere = 1,
                     NumInfectedCritical = 5;
    using InfState                       = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms,
                                            NumInfectedSevere, NumInfectedCritical, 1, 1>;
    using Model    = mio::lsecir::Model<LctState>;
    Model model;

If we do not want to use the default parameters, they can be set as follows.

.. code-block:: cpp

    model.parameters.get<mio::lsecir::TimeExposed>()[0]            = 3.2;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[0] = 2.;
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[0]   = 5.8;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()[0]     = 9.5;
    model.parameters.get<mio::lsecir::TimeInfectedCritical>()[0]   = 7.1;

    model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[0] = 0.05;

    model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[0] = 0.7;
    model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[0] = 0.25;
    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[0] = 0.09;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[0]      = 0.2;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()[0]              = 0.25;
    model.parameters.get<mio::lsecir::DeathsPerCritical>()[0]              = 0.3;

Here, we set the contact matrix used in the simulation. One can define multiple matrices for different locations. The size of each of these matrices is defined by the number of age groups. 
Below, we use only one contact matrix for one location. In our example, we only consider one age group and set the contact rate to 10. 

.. code-block:: cpp

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));

To simulate the implementation of nonpharmaceutical interventions, we add dampings to the contact rate. Here, we apply a damping of :math:`0.7` after :math:`5` days, meaning that the contact rate is reduced to :math:`30%` of the initial value. 

.. code-block:: cpp

    contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.));

For the simulation, we need initial values for all (sub)compartments. If we do not set the initial values manually, these are internally set to :math:`0`.

We start with constructing a vector ``initial_populations`` that we will pass on to the model. It contains vectors for each compartment, that contains a vector with initial values for the respective subcompartments. 
    
.. code-block:: cpp

        std::vector<std::vector<ScalarType>> initial_populations = {{750}, {30, 20},          {20, 10, 10}, {50},
                                                                    {50},  {10, 10, 5, 3, 2}, {20},         {10}};

We assert that vector has the correct size by checking that the number of ``InfectionStates`` and the numbers of subcompartments are correct.

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

We can simulate using the defined model from :math:`t_0` to :math:`t_{\max}` with initial step size :math:`dt` as follows:

.. code-block:: cpp

    ScalarType t0 = 0;
    ScalarType tmax = 10;
    ScalarType dt = 0.5;
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(t0, tmax, dt, model);

The simulation result is divided by subcompartments. We can call the function ``calculate_compartments()`` to get a result according to the ``InfectionStates`` .

.. code-block:: cpp

    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_compartments(result);

We can interpolate the simulation results to a ``TimeSeries`` containing only full days and print the results to the terminal. 

.. code-block:: cpp

    auto interpolated_results = mio::interpolate_simulation_result(population_no_subcompartments);
    interpolated_results.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 12, 4);


Remarks
~~~~~~~~

Above, we have defined the vector of initial values ``initial_populations`` directly. There also exists a function, that computes an intial value vector for the compartments based on a ``TimeSeries`` with flows that are given for a big enough time window before the simulation start. We will demonstarte this below. 
Here, we assume that a model was already constructedas above. 

We start with defining the vectors ``total_population``, ``deaths`` and ``total_confirmed_cases`` that contain the respective values per age group.

.. code-block:: cpp

        Eigen::VectorX<ScalarType> total_population      = Eigen::VectorX<ScalarType>::Constant(1, 1000000.);
        Eigen::VectorX<ScalarType> deaths                = Eigen::VectorX<ScalarType>::Constant(1, 10.);
        Eigen::VectorX<ScalarType> total_confirmed_cases = Eigen::VectorX<ScalarType>::Constant(1, 16000.);



Now, we will define a time series containing flows for some time before the simulation start that will later be used to compute the initial values for the compartments. 

We start by defining the time step size :math:`dt` that determines the distance between the time points that will be added to the time series.  

.. code-block:: cpp

        ScalarType dt = 0.001;

We proceed by creating a time series ``flows`` that contains a vector with the size of the number of transitions that the model allows. 

.. code-block:: cpp
    
        int num_transitions = (int)mio::lsecir::InfectionTransition::Count;
        mio::TimeSeries<ScalarType> flows(num_transitions);

Here, we define the vector that will be added to the time series for each time point. 

.. code-block:: cpp

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


We add the first time point at :math:`-110` and add time points until time :math:`0` where the time step size :math:`dt`determines the distance between the time points. 

.. code-block:: cpp

        flows.add_time_point(-110, vec_flows);
        while (flows.get_last_time() < -dt / 2) {
            flows.add_time_point(flows.get_last_time() + dt, vec_flows);
        }

Now, we can construct an object of type ``Initializer`` where the computations for the initial value vector will be performed.

.. code-block:: cpp

        mio::lsecir::Initializer<Model> initializer(std::move(flows), model);

Finally, we can compute the initialization vector. This is based on the knowledge of the flows as well as the Erlang-distributed stay times in the respective compartments. For further details, see the documentation of the function.

.. code-block:: cpp

        auto status = initializer.compute_initialization_vector(total_population, deaths, total_confirmed_cases);
    


List of models
-----------------------

.. toctree::
    :titlesonly:
    
    models/lsecir
