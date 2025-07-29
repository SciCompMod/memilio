IDE models
==========

MEmilio implements two models based on integro-differential equations (IDEs) with different infection states. IDE-based 
models are a generalization of ODE-based models. Whereas ODE-based models assume an exponential distribution regarding 
the time spent in an infection state, IDE-based models allow for arbitrary distributions.

In the following, we present the general structure of the IDE-SECIR model. 
The IDE-SEIR model follows a different structure, an introduction and a detailed example are provided in a separate 
section below. 

An overview of nonstandard but often used data types can be found under :doc:`data_types`.


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
compartment or the contact rates between different sociodemographic groups. Most model parameters are constants that describe 
pathogen-specific characteristics (possibly resolved by sociodemographic groups) and are represented by a vector with a
value for each group. To model different contact rates between different sociodemographic groups, we
use a parameter denoted **ContactPatterns** of type **UncertainContactMatrix**. The **UncertainContactMatrix** contains an
arbitrary large set of contact matrices and which can represent the different contact locations in the model like 
schools, workplaces, or homes. The matrices can be loaded or stored in the particular example.

In the **ContactPatterns**, each matrix element stores baseline contact rates :math:`c_{i,j}` between sociodemographic 
group :math:`i` and group :math:`j`. The dimension of the matrix is automatically defined by the model initiation and it is reduced 
to one value if no stratification is used. The values can be adjusted during the simulation, e.g., through implementing 
nonpharmaceutical interventions, see the section on :ref:`Nonpharmaceutical Interventions IDE`. 

An important feature of our IDE-based model is that we can choose the transition distributions in a flexible way. The 
default distribution is a smoother cosine function as it provides good testing qualities. For more realistic simulations, 
MEmilio provides the possibility to use exponential, gamma or lognormal distributions within the model.
Practically, one first needs to create an object of a class that is derived from the class ``StateAgeFunction``, 
e.g. ``SmootherCosine``. Any class that is derived from ``StateAgeFunction`` can be inserted into a 
``StateAgeFunctionWrapper`` object that is then passed to the model.

Parameters can get accessed via ``model.parameters.get<Param<ScalarType>>()`` and set via either 
``model.parameters.get<Param<ScalarType>>() = value`` or ``model.parameters.set<Param<ScalarType>>(value)``.


Initial conditions
------------------

The initial conditions consist of a **TimeSeries** that contains the flows per time step for some time interval before 
the simulation start. The length of this time interval depends on the chosen transition distributions and takes into 
account how long there is a relevant fraction of individuals remaining in the compartments. For more information, see 
the documentation of **StateAgeFunction**. Note that the last time point of the initial **TimeSeries** determines the 
start time of the simulation. 


.. _Nonpharmaceutical Interventions IDE:
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
--------------

.. toctree::
    :titlesonly:
    
    models/isecir



IDE-SEIR model
--------------

The IDE-SEIR roughly follows the same structure but differs in a few points:

- We do not consider infection state transitions but only infection states. 
- The initial conditions are defined by the compartment sizes at the simulation start :math:`t_0`.

Note that the IDE-SEIR model is solved with a different numerical scheme than the IDE-SECIR model. For further information, 
see :doc:`IDE-SEIR <./models/iseir>`.

List of models
--------------

.. toctree::
    :titlesonly:
    
    models/iseir




