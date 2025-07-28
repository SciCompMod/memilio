IDE-based SEIR-type model
=========================

The IDE-based SEIR type module models and simulates an epidemic using integro-differential equations. The model is 
particularly suited for pathogens with a latent infection state. The model assumes perfect immunity after recovery and is 
thus only suited for epidemic use cases. In the following, we present the model in detail. Note that this model is 
solved differently tha the IDE-SECIR model and the structure presented in :doc:`IDE models <../ide>` applies only partially.

The simulation runs in discrete time steps using a trapezoidal rule. The model and the numerical scheme is based on the paper 
Keimer A, Pflug L. (2020). *Modeling infectious diseases using integro-differential equations: Optimal
control strategies for policy decisions and Applications in COVID-19*. Preprint. 
`http://dx.doi.org/10.13140/RG.2.2.10845.44000 <http://dx.doi.org/10.13140/RG.2.2.10845.44000>`_. 

For a detailed description and application of the model, see:
Plötzke L. (2020). *Modellierung epidemischer Infektionskrankheiten auf der Basis von gewöhnlichen und Integro-Differentialgleichungen*. 
Bachelor's thesis, University of Cologne. `https://elib.dlr.de/143504/ <https://elib.dlr.de/143504/>`_. 


Infection States
----------------

The model contains the following list of **InfectionState**\s:

.. code-block:: RST

    `Susceptible`
    `Exposed`
    `Infected`
    `Recovered`


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
     - ``ContactFrequency``
     - Average number of contacts of an individual per day.
   * - :math:`\rho`
     - ``TransmissionRisk``
     - Transmission risk for people located in the susceptible compartments.
   * - :math:`T_{E}`
     - ``LatencyTime``
     - Time in days an individual stays in the Exposed compartment.
   * - :math:`T_{I}`
     - ``InfectiousTime``
     - Time in days an individual stays in the Infected compartment.


Initial conditions
------------------

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

.. _Nonpharmaceutical Interventions:
Nonpharmaceutical Interventions
-------------------------------

In the SEIR model, nonpharmaceutical interventions (NPIs) are implemented through dampings in the contact matrix. 
These dampings reduce the contact rates between different groups to simulate interventions.

Basic dampings can be added to the contact matrix as follows:

.. code-block:: cpp

    // Create a contact matrix with constant contact rates between all groups.
    ScalarType cont_freq = 10.;
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    
    // Add a damping that reduces contacts by 70% starting at day 30.
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));


Simulation
----------

After defining :math:`t_{\max}`, we can simulate, which means that we calculate the value for the compartment :math:`S`.

.. code-block:: cpp

    int tmax  = 15;
    model.simulate(tmax);

The values of the remaining compartments :math:`E`, :math:`I` and :math:`R` are calculated using the parameters 
``LatencyTime`` and ``InfectiousTime`` and obtain a time series containing the values of all compartments. 

.. code-block:: cpp

    auto result = model.calculate_EIR();


Output
------

The output of the simulationis a `TimeSeries` object containing the size of the compartments at all time 
points. You can access the results as follows:

.. code-block:: cpp

    // Get the number of time points.
    auto num_points = static_cast<size_t>(result.get_num_time_points());
    
    // Access data at a specific time point.
    Eigen::VectorX value_at_time_i = result.get_value(i);
    ScalarType time_i = result.get_time(i);
    
    // Access the last time point.
    Eigen::VectorX last_value = result.get_last_value();
    ScalarType last_time = result.get_last_time();

The order of the compartments follows the definition in the `InfectionState` enum.

You can print the simulation results as a formatted table:

.. code-block:: cpp

    // Print results to console with default formatting.
    result.print_table();
    
    // Print with custom column labels.
    std::vector<std::string> labels = {"S", "E", "I", "R"};
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

- `examples/ide_seir.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ide_seir.cpp>`_


Overview of the ``iseir`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::iseir