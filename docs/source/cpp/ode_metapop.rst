ODE-based metapopulation model
==============================

Introduction
----------------

This metapopulation model incorporates the effects of spatial heterogeneity and population dynamics into a system of ordinary differential equations (ODEs). A population is divided into several subpopulations, each representing spatially seperated regions which interact on some level. Each subpopulation is further divided into epidemiological compartments and eventually into age groups.
Commuting between regions is governed by a commuting matrix, which describes the fraction of individuals that commute from one region to another. In comparison to the :doc:`Graph-based metapopulation model<graph_metapop>`, commuting is not performed explicitly, i.e., individuals are not exchanged between subpopulations, but rather their theoretical impact to transmission dynamics in other locations is considered, leading to a large system of ODEs which can be solved using standard numerical integration methods.

Simulation
----------

As the model is given as a large system of ODEs, the simulation can be performed by solving the system using a standard ODE solver. The model is build on the same setup as the simpler ODE-models, so we refer to :doc:`ode` for more details.


How to: Set up and run the ODE-based metapopulation model
---------------------------------------------------------

To set up a simulation of the ODE metapopulation model, you need to initialize the model with the desired number of regions and age groups, e.g., 3 regions and 1 age group:
.. code-block:: cpp

    mio::oseirmetapop::Model<double> model(3, 1)

Set a population with the number of individuals in each region, age group and epidemiological compartment, e.g.:
.. code-block:: cpp

    model.populations[{mio::oseirmetapop::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Susceptible}] = 900;
    model.populations[{mio::oseirmetapop::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Exposed}]     = 100;
    model.populations[{mio::oseirmetapop::Region(1), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Susceptible}] = 1000;
    model.populations[{mio::oseirmetapop::Region(2), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Susceptible}] = 1000;

and the epidemiological parameters, e.g.:
.. code-block:: cpp

    model.parameters.template get<ContactPatterns<>>().get_cont_freq_mat()[0].get_baseline().setConstant(2.7);

    model.parameters.set<mio::oseirmetapop::TimeExposed<>>(3.335);
    model.parameters.set<mio::oseirmetapop::TimeInfected<>>(8.097612257);
    model.parameters.set<mio::oseirmetapop::TransmissionProbabilityOnContact<>>(0.07333);

Construct an ``Eigen::MatrixXd`` of size :math:`n_{regions} \times n_{regions}` which describes the fraction of individuals commuting from one region to another. The matrix should satify the sum of each row equal to 1.0, e.g.:
.. code-block:: cpp

    Eigen::MatrixXd mobility_data_commuter(3, 3);
    mobility_data_commuter << 0.4, 0.3, 0.3, 0.2, 0.7, 0.1, 0.4, 0.1, 0.5;

Set the commuting strengths matrix via the ``set_commuting_strengths`` method to ensure that the population after commuting is correctly updated:
.. code-block:: cpp
    
    model.set_commuting_strengths(mobility_data_commuter);

If no mobility data is provided, you should still run ``model.set_commuting_strengths()`` which sets the matrix to an identity matrix, reflecting no mobility between regions. Without setting this, you will encounter errors during the simulation.

Finally, to run the simulation from `t0` to `tmax` with a time step of `dt`, use the following command:
.. code-block:: cpp

    simulate(t0, tmax, dt, model);