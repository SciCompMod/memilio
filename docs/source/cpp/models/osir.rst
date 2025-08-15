ODE-based SIR-type model
=========================

The ODE-based SIR-type module models and simulates the epidemic using an ODE-based SIR-type model approach. The model is
particularly suited for simple simulations of infectious diseases in a population and getting started with the MEmilio 
library. The model assumes perfect immunity after recovery and is thus only suited for epidemic use cases. 
In the following, we present the model in detail.

The infection states and the transitions are visualized in the following image.

.. image:: https://github.com/SciCompMod/memilio/assets/69154294/01c9a2ae-2f5c-4bad-b7f0-34de651f2c73
   :alt: SIR_model


Infection States 
----------------

The model contains the following list of **InfectionState**\s:

.. code-block:: RST

   `Susceptible`
   `Infected`
   `Recovered`


Infection State transitions
---------------------------

The ODE-SIR model is implemented as a **CompartmentalModel**. In each time step, the model computes the aggregated
compartment values.


Sociodemographic Stratification
-------------------------------

In the ODE-SIR model, the population can be stratified by one sociodemographic dimension. This dimension is denoted 
**AgeGroup** but can also be used for other interpretations. For stratification with two or more dimensions, see 
:doc:`Model Creation <../ode_creation>`.

The number of age groups is specified in the model constructor and the model can be initialized with:

.. code-block:: cpp

  mio::oseir::Model model(nb_groups)


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
     - Transmission risk for people located in the Susceptible compartment.
   * - :math:`N`
     - ``populations.get_total()``
     - Total population.
   * - :math:`T_{I}`
     - ``TimeInfected``
     - Time in days an individual stays in the Infected compartment.


Initial conditions
------------------

The initial conditions of the model are defined by the class **Populations** which defines the number of individuals in
each sociodemographic group and **InfectionState**. Before running a simulation, you need to set the initial values for
each compartment:

.. code-block:: cpp

   // Set total population size
   model.populations.set_total(nb_total_t0);

   // Set values for each InfectionState in the specific age group
   model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Infected}] = nb_inf_t0;
   model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Recovered}] = nb_rec_t0;

   // Set the susceptible population as difference to ensure correct total population
   model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osir::InfectionState::Susceptible}, nb_total_t0);

For age-resolved simulations, you need to set the initial conditions for each age group. Additionally, you can use 
``set_difference_from_group_total`` to set the susceptible compartment as the difference between the total group size 
and all other compartments.

.. code-block:: cpp

   for (auto i = mio::AgeGroup(0); i < nb_groups; i++){
      model.populations[{i, mio::osir::InfectionState::Infected}] = 1/nb_groups * nb_inf_t0;
      model.populations[{i, mio::osir::InfectionState::Recovered}] = 1/nb_groups * nb_rec_t0;

      model.populations.set_difference_from_group_total<mio::AgeGroup>(
          {i, mio::osir::InfectionState::Susceptible}, 1/nb_groups * nb_total_t0);
   }


Nonpharmaceutical Interventions
-------------------------------

In the ODE-SIR model, nonpharmaceutical interventions (NPIs) are implemented through dampings in the contact matrix.
These dampings reduce the contact rates between different sociodemographic groups to simulate interventions.

Basic dampings can be added to the ContactPatterns as follows:

.. code-block:: cpp

   contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(nb_groups, nb_groups, 1/nb_groups * cont_freq));

   // Add a uniform damping across all age groups
   contact_matrix.add_damping(Eigen::MatrixXd::Constant(nb_groups, nb_groups, 0.7), mio::SimulationTime(30.));


Simulation
----------

The ODE-SIR model can be simulated using the **simulate** function, which tracks the compartment sizes over
time:

.. code-block:: cpp

   double t0 = 0; // Start time
   double tmax = 50; // End time
   double dt = 0.1; // Time step

   // Run a standard simulation
   mio::TimeSeries<double> result = mio::simulate(t0, tmax, dt, model);

You can also specify a custom integrator:

.. code-block:: cpp

   auto integrator = std::make_shared<mio::RKIntegratorCore>();
   integrator->set_dt_min(0.3);
   integrator->set_dt_max(1.0);
   integrator->set_rel_tolerance(1e-4);
   integrator->set_abs_tolerance(1e-1);

   mio::TimeSeries<double> result = mio::simulate(t0, tmax, dt, model, integrator);


Output
------

The output of the simulation is a ``TimeSeries`` object containing the sizes of each compartment at each time point. For
a basic simulation, you can access the results as follows:

.. code-block:: cpp

   // Get the number of time points
   auto num_points = static_cast<size_t>(result.get_num_time_points());

   // Access data at specific time point 
   Eigen::VectorXd value_at_time_i = result.get_value(i);
   double time_i = result.get_time(i);

   // Access the last time point
   Eigen::VectorXd last_value = result.get_last_value();
   double last_time = result.get_last_time();

You can print the simulation results as a formatted table:

.. code-block:: cpp

   // Print results to console with default formatting
   result.print_table();

   // Print with custom column labels
   std::vector<std::string> labels = {"S", "I", "R"};
   result.print_table(labels);

Additionally, you can export the results to a CSV file:

.. code-block:: cpp

   // Export results to CSV with default settings
   result.export_csv("simulation_results.csv");


Visualization
-------------

To visualize the results of a simulation, you can use the Python package :doc:`m-plot <../../python/m-plot>` 
and its documentation.


Examples
--------

An example can be found at `examples/ode_sir.cpp <https://github.com/SciCompMod/memilio/tree/main/cpp/examples/ode_sir.cpp>`_.


Overview of the ``osir`` namespace:
-----------------------------------------

.. doxygennamespace:: mio::osir
