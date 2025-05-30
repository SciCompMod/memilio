Generalized Linear Chain Trick model
===========

MEmilio implements a SECIR-type model utilizing the  Generalized Linear Chain Trick (GLCT). This is a generalization of 
the LCT model and allows for phase-type distributed stay times in the compartments. Phase-type distributions are dense 
in the field of all positive-valued distributions. Therefore, for any positive-valued distribution, a phase-type 
distribution of arbitrary precision can be identified. Note that the resulting system can still be described by ODEs. 

In the following, we present the general structure of the GLCT model. The particular model documentation with examples 
is linked at the bottom of this page.

Infection states
----------------

The model contains a list of **InfectionState**\s that define particular features of the subpopulations in the particular state.

.. code-block:: RST

    `State1`
    `State2`
    `...`

To make use of the GLCT, we additionally need to define the numbers of subcompartments for each **InfectionState**.

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

In the following, we give detailed explanations of the GLCT-SECIR model.

Introduction
-------------

This model is based on the Generalized Linear Chain Trick (GLCT). 

The GLCT provides the option to use phase-type distributed stay times in the compartments through the use of subcompartments. The Generalized Linear Chain Trick is an extension of the Linear Chain Trick (as the name already suggests). Phase-type distributions are dense in the field of all positive-valued distributions. Therefore, for any positive-valued distribution, a phase-type distribution of arbitrary precision can be identified.
The normal ODE models have (possibly unrealistic) exponentially distributed stay times.
The GLCT model can still be described by an ordinary differential equation system.

For the concept see 

- Hurtado PJ and Kirosingh AS (2019) Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models. Journal of Mathematical Biology. https://doi.org/10.1007/s00285-019-01412-w
- Hurtado PF and Richards C (2021) Building mean field ODE models using the generalized linear chain trick & Markov chain theory. Journal of Biological Dynamics. https://doi.org/10.1080/17513758.2021.1912418

Here, the eight compartments 

- `Susceptible` (:math:`S`), may become Exposed at any time
- `Exposed` (:math:`E`), becomes InfectedNoSymptoms after some time
- `InfectedNoSymptoms` (:math:`I_{NS}`), becomes InfectedSymptoms or Recovered after some time
- `InfectedSymptoms` (:math:`I_{Sy}`), becomes InfectedSevere or Recovered after some time
- `InfectedSevere` (:math:`I_{Sev}`), becomes InfectedCritical or Recovered after some time
- `InfectedCritical` (:math:`I_{Cr}`), becomes Recovered or Dead after some time
- `Recovered` (:math:`R`)
- `Dead` (:math:`D`)

are used to simulate the spread of the disease. 
It is possible to include phase-type distributed stay times for the five compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.

Simulation
-----------

The simulation runs in discrete time steps. Different ODE solvers are available, some of them use an adaptive time step size.

For a detailed description of the concept, see:

- Hurtado PJ und Kirosingh AS (2019) Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models. Journal of Mathematical Biology. https://doi.org/10.1007/s00285-019-01412-w
- Hurtado PJ und Richards C (2021) Building mean field ODE models using the generalized linear chain trick & Markov chain theory. Journal of Mathematical Biology. https://doi.org/10.1080/17513758.2021.1912418

How to: Set up and run a simulation of the GLCT-SECIR model
------------------------------------------------------------

In the following, we will demonstrate how to run a simulation using the GLCT-SECIR model. The here defined example recreates the lct_secir.cpp example using the GLCT model. This means, that we use the corresponding initial numbers, parameters and numbers of subcompartments that are equivalent to the choices in the LCT example.

We start by defining the number of subcompartments and constructing the model. We can choose the number of subcompartments individually for the compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.
Note that in the GLCT model, we define two strains for the compartments `InfectedNoSymptoms`, `InfectedSymptoms`, `InfectedSevere` and `InfectedCritical` as individuals in these compartments can either transition to an infection state corresponding to a more severe disease state or recover. This is why we define the model with twice the number of subcompartments compared to the LCT-SECIR model for these infection states. 

.. code-block:: cpp

    constexpr size_t NumExposed = 2, NumInfectedNoSymptoms = 6, NumInfectedSymptoms = 2, NumInfectedSevere = 2,
                    NumInfectedCritical = 10;
    using Model    = mio::glsecir::Model<NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms, NumInfectedSevere,
                                    NumInfectedCritical>;
    using LctState = Model::LctState;
    using InfectionState = LctState::InfectionState;

    Model model;

We continue by defining some epidemiological parameters needed throughout the model definition and initialization.

.. code-block:: cpp

    const ScalarType timeExposed                    = 3.2;
    const ScalarType timeInfectedNoSymptoms         = 2.;
    const ScalarType timeInfectedSymptoms           = 5.8;
    const ScalarType timeInfectedSevere             = 9.5;
    const ScalarType timeInfectedCritical           = 7.1;
    const ScalarType recoveredPerInfectedNoSymptoms = 0.09;
    const ScalarType severePerInfectedSymptoms      = 0.2;
    const ScalarType criticalPerSevere              = 0.25;
    const ScalarType deathsPerCritical              = 0.3;

Now, we define the initial values with the distribution of the population into subcompartments. Note that this method of defining the initial values using a vector of vectors is not necessary, but should remind you how the entries of the initial value vector relate to the defined template parameters of the model or the number of subcompartments. It is also possible to define the initial values directly.

In this example, we want to initialize the GLCT model so that it corresponds to the example given for the LCT model. This is why we take the initial population from the LCT example and split it into two strains according to the respective transition probabilities for the the compartments InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.
Note that in the case of InfectedNoSymptoms, the first three subcompartments correspond to the first strain, i.e. the individuals that will transition to InfectedSymptoms afterwards, and the other three subcompartments correspond to the second strain, i.e. the individuals that willl recover. For the other compartments for which we defined two strains in the model, this is done analogously. 

.. code-block:: cpp

    std::vector<std::vector<ScalarType>> initial_populations = {
        {750}, // Susceptible
        {30, 20}, // Exposed
        {20 * (1 - recoveredPerInfectedNoSymptoms), 10 * (1 - recoveredPerInfectedNoSymptoms), // InfectedNoSymptoms
        10 * (1 - recoveredPerInfectedNoSymptoms), 20 * recoveredPerInfectedNoSymptoms,
        10 * recoveredPerInfectedNoSymptoms, 10 * recoveredPerInfectedNoSymptoms},
        {50 * severePerInfectedSymptoms, 50 * (1 - severePerInfectedSymptoms)}, // InfectedSymptoms
        {50 * criticalPerSevere, 50 * (1 - criticalPerSevere)}, // InfectedSevere
        {10 * deathsPerCritical, 10 * deathsPerCritical, 5 * deathsPerCritical, 3 * deathsPerCritical, // InfectedCritical
        2 * deathsPerCritical, 10 * (1 - deathsPerCritical), 10 * (1 - deathsPerCritical), 5 * (1 - deathsPerCritical),
        3 * (1 - deathsPerCritical), 2 * (1 - deathsPerCritical)},
        {20}, // Recovered
        {10}}; // Dead

Here, we assert that ``initial_populations`` has the right shape.

.. code-block:: cpp

    if (initial_populations.size() != (size_t)InfectionState::Count) {
        mio::log_error("The number of vectors in initial_populations does not match the number of InfectionStates.");
        return 1;
    }
    if ((initial_populations[(size_t)InfectionState::Susceptible].size() !=
        LctState::get_num_subcompartments<InfectionState::Susceptible>()) ||
        (initial_populations[(size_t)InfectionState::Exposed].size() != NumExposed) ||
        (initial_populations[(size_t)InfectionState::InfectedNoSymptoms].size() != NumInfectedNoSymptoms) ||
        (initial_populations[(size_t)InfectionState::InfectedSymptoms].size() != NumInfectedSymptoms) ||
        (initial_populations[(size_t)InfectionState::InfectedSevere].size() != NumInfectedSevere) ||
        (initial_populations[(size_t)InfectionState::InfectedCritical].size() != NumInfectedCritical) ||
        (initial_populations[(size_t)InfectionState::Recovered].size() !=
        LctState::get_num_subcompartments<InfectionState::Recovered>()) ||
        (initial_populations[(size_t)InfectionState::Dead].size() !=
        LctState::get_num_subcompartments<InfectionState::Dead>())) {
        mio::log_error("The length of at least one vector in initial_populations does not match the related number of "
                    "subcompartments.");
        return 1;
    }

Finally, we transfer the initial values in ``initial_populations`` to the model.

.. code-block:: cpp

    std::vector<ScalarType> flat_initial_populations;
    for (auto&& vec : initial_populations) {
        flat_initial_populations.insert(flat_initial_populations.end(), vec.begin(), vec.end());
    }
    for (size_t i = 0; i < LctState::Count; i++) {
        model.populations[mio::Index<LctState>(i)] = flat_initial_populations[i];
    }

Since we want to recreate the LCT model as defined in the corresponding example, we want to set the parameters determining the transition distributions such that we obtain Erlang distributions. 

We will explain how to do this for the different compartments below. In general, we need to define a vector ``StartingProbabilities(...)`` and a matrix ``TransitionMatrix(...z)To(...*)`` for all compartments with subcompartments, i.e. `Exposed`, `InfectedNoSymptoms`, `InfectedSymptoms`, `InfectedSevere` and `InfectedCritical`.

The size of the vector ``StartingProbabilities`` is the number of subcompartments and contains the initial probability of starting in any of the subcompartments of the respective compartment. The entries should sum up to 1.

The matrix ``TransitionMatrix(...z)To(...*)`` describes the transitions in between of the subcompartments of compartment `z` to the compartment `*`.


We start with the `Exposed` compartment. The folllowing definitions of the starting probabilities and the transition matrix lead to an Erlang-distributed latent stage.

The get_default of the ``StartingProbabilities(...)`` returns the first unit vector of the defined size. It is necessary to set it although the default method is used to define the length of the vector.

.. code-block:: cpp

    model.parameters.get<mio::glsecir::StartingProbabilitiesExposed>() =
        mio::glsecir::StartingProbabilitiesExposed().get_default(
            LctState::get_num_subcompartments<InfectionState::Exposed>());

The get_default function returns the ``TransitionMatrix`` that is required to have an Erlang-distributed stay time with an average of timeExposed.

.. code-block:: cpp

    model.parameters.get<mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms>() =
        mio::glsecir::TransitionMatrixExposedToInfectedNoSymptoms().get_default(
            LctState::get_num_subcompartments<InfectionState::Exposed>(), timeExposed);
    

We continue with the compartment `InfectedNoSymptoms`. For InfectedNoSymptoms, two strains have to be defined, one for the transition `InfectedNoSymptomsToInfectedSymptoms` and one for the transition `InfectedNoSymptomsToRecovered`.
The strains have a length of ``NumInfectedNoSymptoms/2`` each as we choose the same number of subcompartments for both strains. Note that the transition probability is included in the vector ``StartingProbabilitiesInfectedNoSymptoms``.

.. code-block:: cpp

    Eigen::VectorX<ScalarType> StartingProbabilitiesInfectedNoSymptoms =
        Eigen::VectorX<ScalarType>::Zero(LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>());
    StartingProbabilitiesInfectedNoSymptoms[0] = 1 - recoveredPerInfectedNoSymptoms;
    StartingProbabilitiesInfectedNoSymptoms[(Eigen::Index)(
        LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2.)] = recoveredPerInfectedNoSymptoms;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedNoSymptoms>() =
        StartingProbabilitiesInfectedNoSymptoms;

Define equal transition matrices for the strains. They follow the same Erlang distribution such that we get the same result as with LCT model that can only consider one strain.

.. code-block:: cpp

    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>() =
        mio::glsecir::TransitionMatrixInfectedNoSymptomsToInfectedSymptoms().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2.),
            timeInfectedNoSymptoms);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedNoSymptomsToRecovered().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedNoSymptoms>() / 2.),
            timeInfectedNoSymptoms);

We proceed analogously for the remaining compartments `InfectedSymptoms`, `InfectedSevere` `InfectedCritical`.

.. code-block:: cpp

    // InfectedSymptoms.
    Eigen::VectorX<ScalarType> StartingProbabilitiesInfectedSymptoms =
        Eigen::VectorX<ScalarType>::Zero(LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>());
    StartingProbabilitiesInfectedSymptoms[0]                                         = severePerInfectedSymptoms;
    StartingProbabilitiesInfectedSymptoms[(Eigen::Index)(
        LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2.)] = 1 - severePerInfectedSymptoms;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSymptoms>() = StartingProbabilitiesInfectedSymptoms;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToInfectedSevere().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2.), timeInfectedSymptoms);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedSymptomsToRecovered().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSymptoms>() / 2.), timeInfectedSymptoms);

    // InfectedSevere.
    Eigen::VectorX<ScalarType> StartingProbabilitiesInfectedSevere =
        Eigen::VectorX<ScalarType>::Zero(LctState::get_num_subcompartments<InfectionState::InfectedSevere>());
    StartingProbabilitiesInfectedSevere[0]                                         = criticalPerSevere;
    StartingProbabilitiesInfectedSevere[(Eigen::Index)(
        LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2.)] = 1 - criticalPerSevere;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedSevere>() = StartingProbabilitiesInfectedSevere;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical>() =
        mio::glsecir::TransitionMatrixInfectedSevereToInfectedCritical().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2.), timeInfectedSevere);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedSevereToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedSevereToRecovered().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedSevere>() / 2.), timeInfectedSevere);

    // InfectedCritical.
    Eigen::VectorX<ScalarType> StartingProbabilitiesInfectedCritical =
        Eigen::VectorX<ScalarType>::Zero(LctState::get_num_subcompartments<InfectionState::InfectedCritical>());
    StartingProbabilitiesInfectedCritical[0]                                         = deathsPerCritical;
    StartingProbabilitiesInfectedCritical[(Eigen::Index)(
        LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.)] = 1 - deathsPerCritical;
    model.parameters.get<mio::glsecir::StartingProbabilitiesInfectedCritical>() = StartingProbabilitiesInfectedCritical;
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToDead>() =
        mio::glsecir::TransitionMatrixInfectedCriticalToDead().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.), timeInfectedCritical);
    model.parameters.get<mio::glsecir::TransitionMatrixInfectedCriticalToRecovered>() =
        mio::glsecir::TransitionMatrixInfectedCriticalToRecovered().get_default(
            (size_t)(LctState::get_num_subcompartments<InfectionState::InfectedCritical>() / 2.), timeInfectedCritical);

We set the remaining parameters of the model ``TransmissionProbabilityOnContact``, ``RelativeTransmissionNoSymptoms`` and ``RiskOfInfectionFromSymptomatic``.

.. code-block:: cpp

    model.parameters.get<mio::glsecir::TransmissionProbabilityOnContact>() = 0.05;
    model.parameters.get<mio::glsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    model.parameters.get<mio::glsecir::RiskOfInfectionFromSymptomatic>() = 0.25;

Here, we set the contact matrix used in the simulation. One can define multiple matrices for different locations. The size of each of these matrices is defined by the number of age groups. 
In our example below we use only one contact matrix for one location. Our model considers one age group and we set the contact rate to 10. 

.. code-block:: cpp

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::glsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));

To simulate the implementation of nonpharmaceutical interventions, we add dampings to the contact rate. Here, we apply a damping of :math:`0.7` after :math:`5`` days, meaning that the contact rate is reduced to 30% of the initial value. 

.. code-block:: cpp

contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.));

We can simulate using the defined model from :math:`t_0` to :math:`t_{\max}` with initial step size :math:`dt_init` as follows:

.. code-block:: cpp

    const ScalarType t0      = 0;
    const ScalarType tmax    = 10;
    const ScalarType dt_init = 10;
        mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(t0, tmax, dt_init, model);
    
The simulation result is divided by subcompartments as in the LCT model. We call the function calculate_compartments to get a result according to the ``InfectionStates``.

.. code-block:: cpp

    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_compartments(result);

We can interpolate the simulation results to a ``TimeSeries`` containing only full days and print the results to the terminal. 

.. code-block:: cpp
    
    auto interpolated_result = mio::interpolate_simulation_result(population_no_subcompartments, 0.1);
    interpolated_result.print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 10, 4);


List of models
-----------------------------

.. toctree::
    :titlesonly:
    
    models/glsecir
