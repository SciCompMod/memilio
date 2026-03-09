GLCT SECIR model
================

The GLCT-SECIR module models and simulates an epidemic using an ODE-based approach while making use of the Generalized Linear Chain Trick to provide the option of phase-type distributed stay times in the compartments through the use of subcompartments. The model is particularly suited for pathogens with pre- or asymptomatic infection states and when severe or critical states are possible. The model assumes perfect immunity after recovery and is thus only suited for epidemic use cases. In the following, we present the model in detail.

The compartment structure with subcompartments is the same as in the LCT-SECIR model. An overview of the model 
architecture can be found in the `LCT model <lsecir>`_. For the GLCT model, some additional transitions are possible and we have more arrows in the model architecture. Below is an example for the Exposed compartment. Note that some indices are omitted (e.g., :math:`n` instead of :math:`n_E`) to keep the picture simple.

.. image:: https://github.com/user-attachments/assets/fc075b7a-6cd2-4e70-bdd0-a2f4b9f2cf53
   :alt: tikzGLCTSECIR

For the concept see:

- Hurtado PJ, Kirosingh AS. (2019). *Generalizations of the ‘Linear Chain Trick’: incorporating more flexible dwell time distributions into mean field ODE models*. Journal of Mathematical Biology. `https://doi.org/10.1007/s00285-019-01412-w <https://doi.org/10.1007/s00285-019-01412-w>`_ 
- Hurtado PJ, Richards C. (2021). *Building mean field ODE models using the generalized linear chain trick & Markov chain theory*. Journal of Mathematical Biology. `https://doi.org/10.1080/17513758.2021.1912418 <https://doi.org/10.1080/17513758.2021.1912418>`_  

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

Parameters
---------------

Below is an overview of the model variables:

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
     - Transmission risk for people located in the Susceptible compartments.
   * - :math:`\xi_{I_{NS}}`
     - ``RelativeTransmissionNoSymptoms``
     - Proportion of nonsymptomatically infected people who are not isolated.
   * - :math:`\xi_{I_{Sy}}`
     - ``RiskOfInfectionFromSymptomatic``
     - Proportion of infected people with symptoms who are not isolated.
   * - :math:`n_{z}`
     - ``Num(...)``
     - Number of subcompartments of compartment :math:`z \in \mathcal{Z}`. (...) refers to the C++-name of :math:`z` as stated above.
   * - :math:`\boldsymbol{\alpha_{z}}`
     - ``StartingProbabilities(...)``
     - Vector of size :math:`n_{z}` with the initial probability of starting in any of the subcompartments of compartment :math:`z \in \mathcal{Z}`. The entries should sum up to 1.
   * - :math:`\mathbf{A_{z}^{*}}`
     - ``TransitionMatrix(...z)To(...*)``
     - Matrix describing the transitions in between the subcompartments of :math:`z \in \mathcal{Z}` that describes the transition to the compartment *.

The model equations are given below. For a simpler description let :math:`\mathcal{Z}=\{E,I_{NS},I_{Sy},I_{Sev},I_{Cr}\}` be the set of the compartments that can be divided into subcompartments.

.. image:: https://github.com/SciCompMod/memilio/assets/70579874/e1da5e1d-e719-4c16-9f14-45374be7c353
   :alt: equations

Note that the bold notation :math:`\mathbf{z}(t)` for :math:`z \in \mathcal{Z}` stands for a vector. If several transitions are possible from a compartment, the vector is split in order to be able to select the stay times until the transitions individually. For example, the order

.. math::

   \mathbf{I_{\text{NS}}}(t) = \begin{bmatrix}
   \mathbf{I_{\text{NS}}^{\text{Sy}}}(t) \\
   \mathbf{I_{\text{NS}}^{\text{R}}}(t)
   \end{bmatrix}

is used. Similar holds true for the other compartments :math:`z \in \mathcal{Z}`.

Implicitly, the matrices :math:`\mathbf{A_{z}^{*}}` for one :math:`z \in \mathcal{Z}` are a block of a matrix :math:`\mathbf{A_{z}}` corresponding to the whole vector :math:`\mathbf{z}(t)`. As we have no transitions in between the strains defined for different transition probabilities, we would have many zeros in the matrix. The matrix can be defined as

.. math::

   \mathbf{A_{z}}=
   \begin{bmatrix}
   \mathbf{A_{z}^{*_1}} &  \mathbf{0} \\
   \mathbf{0} &  \mathbf{A_{z}^{*_2}}
   \end{bmatrix},

where :math:`{*}_{1}` is the compartment of the first transition, e.g., :math:`I_{\text{Sy}}` for :math:`z=I_{\text{NS}}` and :math:`*_{2}` the compartment of the second possible transition, e.g., :math:`R`. Therefore, we just store the non-zero blocks of the matrix. Using these parameters, the phase-type distribution that defines the stay time in compartment :math:`z \in \mathcal{Z}` has the probability density function

.. math::

   f(x)=\boldsymbol{\alpha_z}^T\, e^{x\,\mathbf{A_z}}\, \Bigl(-\mathbf{A_z}\,\boldsymbol{\Bbb{1}}\Bigr)
   \quad \text{for } x\in\mathbb{R}^{+}

and the cumulative distribution function

.. math::

   F(x)=1-\boldsymbol{\alpha_z}^T\, e^{x\,\mathbf{A_z}}\, \boldsymbol{\Bbb{1}},

where

.. math::

   e^{x\,\mathbf{A_z}}=\sum_{j=0}^{\infty}\frac{\bigl(x\,\mathbf{A_z}\bigr)^j}{j!}

is the matrix exponential and :math:`\boldsymbol{\Bbb{1}}` is the vector containing ones of the matching size. Therefore, by changing the vector :math:`\boldsymbol{\alpha_z}` and the matrices :math:`\mathbf{A_{z}^{*}}`, one can choose the stay time distribution appropriately.

It is important that the sizes of the vectors and matrices match each other and satisfy some other conditions that are checked before a simulation.


Initial conditions
------------------

We start by defining the number of subcompartments and constructing the model with it. We can choose the number of subcompartments individually for the compartments Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and InfectedCritical.
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

Now, we define the initial values with the distribution of the population into subcompartments. Note that this method of defining the initial values using a vector of vectors is not necessary, but should show how the entries of the initial value vector relate to the defined template parameters of the model or the number of subcompartments. It is also possible to define the initial values directly.

In this example, we initialize the GLCT model so that it corresponds to the example given for the LCT model. 
For that, we take the initial population from the LCT example and split it into two strains according to the 
respective transition probabilities for the compartments InfectedNoSymptoms, InfectedSymptoms, InfectedSevere and 
InfectedCritical.
In the case of InfectedNoSymptoms, the first three subcompartments correspond to the first strain, i.e. the 
individuals that will transition to InfectedSymptoms afterwards, and the other three subcompartments correspond to the second strain, i.e. the individuals that willl recover. For the other compartments for which we defined two strains in the model, this is done analogously. 

We continue by defining some epidemiological parameters needed throughout the model definition and initialization.

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

Below, we assert that ``initial_populations`` has the right shape.

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


Since we want to recreate the LCT model as defined in the corresponding example, we set the parameters determining the transition distributions such that we obtain Erlang distributions. 

We will explain how to do this for the different compartments in the following. In general, we need to define a vector ``StartingProbabilities(...)`` and a matrix ``TransitionMatrix(...z)To(...*)`` for all compartments with subcompartments, i.e. `Exposed`, `InfectedNoSymptoms`, `InfectedSymptoms`, `InfectedSevere` and `InfectedCritical`.

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

Equal transition matrices for the strains have to be defined. They follow the same Erlang distribution such that we get the same result as with the LCT model that can only consider one strain.

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


.. _Nonpharmaceutical Interventions:
Nonpharmaceutical Interventions
-------------------------------

In the GLCT-SECIR model, nonpharmaceutical interventions (NPIs) are implemented through dampings in the contact matrix. 
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



For more complex scenarios, such as real-world venue closures or lockdown modeling, you can implement detailed NPIs with location-specific dampings. The GLCT-SECIR model supports contact matrices for different locations (e.g., home, school, work, other) and can apply different dampings to each location.

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

Simulating the model from :math:`t_0` to :math:`t_{\max}` with initial step size :math:`dt_init` id done as follows:

.. code-block:: cpp

    const ScalarType t0      = 0;
    const ScalarType tmax    = 10;
    const ScalarType dt_init = 10;
        mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(t0, tmax, dt_init, model);

You can also specify a custom integrator:

.. code-block:: cpp

    auto integrator = std::make_unique<mio::RKIntegratorCore>();
    integrator->set_dt_min(0.3);
    integrator->set_dt_max(1.0);
    integrator->set_rel_tolerance(1e-4);
    integrator->set_abs_tolerance(1e-1);
    
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(t0, tmax, dt, model, std::move(integrator));

Output
------

The simulation result is divided by subcompartments. The function ``calculate_compartments()`` aggregates the subcompartments by `InfectionState`\s .

.. code-block:: cpp

    mio::TimeSeries<ScalarType> population_no_subcompartments = model.calculate_compartments(result);

You can access the data in the `mio::TimeSeries` object as follows:

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

To visualize the results of a simulation, you can use the Python package :doc:`m-plot <../../python/m-plot>` and its documentation.

    
Examples
--------

An example can be found at:

- `examples/glct_secir.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/glct_secir.cpp>`_ 


Overview of the ``glsecir`` namespace:
--------------------------------------

.. doxygennamespace:: mio::glsecir
