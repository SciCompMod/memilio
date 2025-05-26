ODE Models
==========

MEmilio implements various models based on ordinary differential equations (ODEs). ODE-based models are a subclass of 
compartmental models in which individuals are grouped into compartments. MEmilio's ODE-based models range from most simple 
SIR structure to complex models with multi-layer and waning immunity. These models can be stratified by age groups or 
other sociodemographic factors.

In the following, we present the general structure for all simple ODE-based models. For generalizations via the Linear 
Chain Trick and Erlang distributed stay times, see :doc:`here <lct>`. For spatially resolved Graph-ODE models, 
see :doc:`here <graph_metapop>`. The particular model documentations with examples are linked at the bottom of this page.

Infection States
----------------

Every model contains a list of **InfectionState**\s that define particular features of the subpopulations in the particular state.

.. code-block:: RST

    `State1`
    `State2`
    `...`


Infection State Transitions
---------------------------

Our ODE-based models are either implemented as **FlowModel** or as **CompartmentalModel**. In a **FlowModel**, flows 
``Flow<State1, State2>`` between **InfectionState**\s **State1** and **State2** are defined. Instead of a standard 
solution to an ODE-based model, this implementation additionally realizes the solution of the transitions or flows 
between states and directly enables users to access new transmissions or hospitalizations at any time point. 
The simpler class **CompartmentalModel** only considers the states of the system and not the flows.


Sociodemographic Stratification
-------------------------------

For most models, the population can also be stratified by one sociodemographic dimension. This dimension is denoted 
**AgeGroup** but can also be used for other interpretations. For stratifications with two or more dimensions, 
see :doc:`Model Creation <ode_creation>`.


Parameters
----------

The parameters of the model are defined as structs and are combined in a class ``ParameterSet<Param1, Param2, ...>``.
We use different types of parameters to represent epidemiological parameters such as the mean stay times in a 
compartment or the contact rates between different age groups. Most model parameters are constants that describe 
pathogen-specific characteristics (possibly resolved by sociodemographic groups) and are represented by a vector with a
value for each sociodemographic group. To model different contact rates between different sociodemographic groups, we
use a parameter of type **ContactMatrix** that can be adjusted during the simulation, e.g., through implementing 
nonpharmaceutical interventions, see the section on :ref:`Nonpharmaceutical Interventions`.
Parameters can get accessed via ``model.parameters.get<Param>()`` and set via ``model.parameters.set<Param>(value)``.


Initial conditions
------------------

The initial conditions of the model are represented by a class **Populations** that gives the number of individuals in 
each sociodemographic group and **InfectionState**. For more details, see :doc:`Model Creation <ode_creation>`. Before 
the simulation, set the initial conditions via ``model.populations[{AgeGroup::Age, InfectionState::State}] = value`` for
each **InfectionState** and sociodemographic group.


.. _Nonpharmaceutical Interventions:
Nonpharmaceutical Interventions
-------------------------------

Contact rates can be adjusted during the simulation to model nonpharmaceutical interventions (NPIs) such as lockdowns, 
school closures, or social distancing. This is done by adding **Damping**\s to the **ContactMatrix** of the model. A 
**Damping** is defined by a time point at which the intervention starts and a matrix of the same size as the 
**ContactMatrix** giving a factor for each combination of sociodemographic groups by which the corresponding contact 
rate is reduced.


Simulation
----------

Once the model is setup, run a simple simulation from time ``t0`` to ``tmax`` with an initial step size ``dt`` using the 
``mio::simulation()`` function. This will run a simulation of type **Simulation** that does not save the flows between 
compartments but only the sizes of each compartment over time. To use the flow information, make sure to use a 
**FlowModel** and run a simulation of type **FlowSimulation** with the ``mio::simulate_flows()`` function.
You can run a simulation using either fixed or adaptive solution schemes with an absolute or relative tolerance. By 
default, the simulation uses an adaptive solution scheme of the boost library and an absolute tolerance of 1e-10 and a 
relative tolerance of 1e-5. For more details on the possible integration schemes, see <numerical integrator stuff>.


Output
------

The output of the **Simulation** is a ``TimeSeries`` containing the sizes of each compartment at each time point. A 
simple table can be printed using the ``print_table()`` function of the ``TimeSeries`` class. The output of the 
**FlowSimulation** additionally contains the flows between compartments at each time point. The compartment sizes can be
printed with ``result[0].print_table()`` and the flows with ``result[1].print_table()``. 
As adaptive step size methods may be used, you can interpolate the results to days or any other series of times points
with ``mio::interpolate_simulation_result()``.


Visualization
-------------

To visualize the results of a simulation, you can use the Python package :doc:`memilio_plot <../python/memilio_plot>`
and its documentation.


List of models
=====================

.. toctree::
    :titlesonly:
    
    models/osir
    models/oseir
    models/oseair
    models/osecir
    models/osecirvvs
    models/osecirts

    