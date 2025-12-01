ODE-based metapopulation model
================================

Introduction
------------

This model incorporates spatial resolution into an ODE-based model by including it into the differential equations. The population under study is divided into multiple regions and potentially into several sociodemographic groups such as age groups. Every subgroup is further divided into the epidemiological compartments Susceptible, Exposed, Infected and Recovered. The mobility is deterministically included into the ODEs, forming a large system of ODEs with an equation for every combination of region, sociodemographic group and epidemiological compartment.


Simulation
----------

As the model is given as a large system of ODEs, the simulation can be performed by solving the system using a standard ODE solver. The model is build on the same setup as the simpler ODE-models, so we refer to :doc:`ode` for more details.


How to: Set up and run the ODE metapopulation model
---------------------------------------------------

To set up a simulation of the ODE metapopulation model, you need to set the parameters and initial populations for each region. To do so, we can proceed similarly as for the simpler ODE models, see :doc:`ode`. Additionally, we need to set up the mobility between the regions via a mobility matrix of size :math:`n\times n`:

.. code-block:: cpp

    Eigen::MatrixXd mobility_data_commuter(n, n);
    mobility_data_commuter << ...; // Fill the matrix with mobility data

Afterwards, it can be set in the model via ``model.set_commuting_strengths(mobility_data_commuter)``. If no mobility data is provided, a default mobility matrix with zero mobility between all regions is used. Note, that it is still necessary to run the ``model.set_commuting_strengths()`` function to initialize the default mobility in the model.

Similar to the simple ODE models, the simple simulation can be run via the ``mio::simulation()`` function. For more details on the simulation, see :doc:`ode`.