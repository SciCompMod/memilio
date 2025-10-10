ODE-based models
================

MEmilio implements various models based on ordinary differential equations (ODEs). ODE-based models are a subclass of 
compartmental models in which individuals are grouped into subpopulations called compartments. MEmilio's ODE-based models range from most simple 
SIR structures to complex models with multiple layers and waning immunity. These models can be stratified by age groups or 
other sociodemographic factors.

In the following, we present the general structure for all simple ODE-based models. For generalizations via the Linear 
Chain Trick and Erlang distributed stay times, see :doc:`here <lct>`. For spatially resolved Graph-ODE models, 
see :doc:`here <graph_metapop>`. The particular model documentations with examples are linked at the bottom of this page.

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

Our ODE-based models are either implemented as **FlowModel** or as **CompartmentalModel**. In a **FlowModel**, flows 
``Flow<State1, State2>`` between **InfectionState**\s **State1** and **State2** are defined. Instead of a standard 
solution of an ODE-based model, this implementation additionally realizes the solution of the transitions or flows 
between states and directly enables users to access new transmissions or hospitalizations at any time point. 
The simpler class **CompartmentalModel** only considers the states of the system and not the flows.


Sociodemographic stratification
-------------------------------

For most ODE-based models, the population can be stratified by one sociodemographic dimension. This dimension is denoted **AgeGroup** but can also be used for other interpretations. For stratifications with two or more dimensions, see :doc:`Model Creation <ode_creation>`.


Parameters
----------

The parameters of the model are defined as structs and are combined in a class ``ParameterSet<Param1, Param2, ...>``.
We use different types of parameters to represent epidemiological parameters such as the mean stay times in a 
compartment or the contact rates between different age groups. Most model parameters are constants that describe 
pathogen-specific characteristics (possibly resolved by sociodemographic groups) and are represented by a vector with a
value for each sociodemographic group. To model different contact rates between different sociodemographic groups, we
use a parameter denoted **ContactPatterns** of type **UncertainContactMatrix**. The **UncertainContactMatrix** contains
a set of contact matrices of arbitrary length which can represent different contact locations in the model like 
schools, workplaces, or homes. The matrices can be loaded or stored in the particular example.
In the **ContactPatterns** parameter, each matrix element stores baseline contact rates :math:`c_{i,j}` between sociodemographic 
group :math:`i` and group :math:`j`. The dimension of the matrix is automatically defined by the model initiation and it is reduced 
to one value if no stratifcation is used. The values can be adjusted during the simulation, e.g., through implementing 
nonpharmaceutical interventions, see the section on :ref:`Nonpharmaceutical Interventions`.
Parameters can get accessed via ``model.parameters.get<Param<double>>()`` and set via either 
``model.parameters.get<Param<double>>() = value`` or ``model.parameters.set<Param<double>>(value)``. 

.. dropdown:: :fa:`gears` Expert's settings

    In the above description, ``double`` has been used for the template parameter ``FP``. ``FP`` could also be replaced by an
    object type for automatic differentiation and, perspectively, should allow the computation in lower precision to 
    reduce the computational runtime.


Initial conditions
------------------

The initial conditions of the model are represented by a class **Populations** that gives the number of individuals in 
each sociodemographic group and **InfectionState**. For more details, see :doc:`Model Creation <ode_creation>`. Before 
the simulation, set the initial conditions via ``model.populations[{AgeGroup::Age, InfectionState::State}] = value`` for
each **InfectionState** and sociodemographic group.


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

Once the model is set up, one can run a simple simulation from time ``t0`` to ``tmax`` with an initial step size ``dt`` using the 
``mio::simulate()`` function. This will run a simulation of type **Simulation** that saves the sizes of each compartment over time. 
To also save the flow information, make sure to use a **FlowModel** and run a simulation of type **FlowSimulation** with the ``mio::simulate_flows()`` function.
You can run a simulation using either fixed or adaptive integration schemes with an absolute or relative tolerance. By 
default, the simulation uses an adaptive integration scheme of the boost library with an absolute tolerance of 1e-10 and a 
relative tolerance of 1e-5.

Additionally, a feedback simulation is available, which allows for dynamic adjustments of the simulation based on its own output. This is particularly useful for modeling behavioral changes in response to the epidemiological situation. The local feedback simulation is based on the work of Dönges et al. (doi.org/10.3389/fphy.2022.842180). It uses the history of ICU occupancy to calculate a "perceived risk", which then translates into a reduction of contact rates. The ``FeedbackSimulation`` class wraps an existing simulation and applies this feedback mechanism at regular intervals. This concept can be extended to a metapopulation model using a ``FeedbackGraphSimulation``, which is described in more detail in the :doc:`Graph-based metapopulation model <graph_metapop>` section.


Output
------

The output of the **Simulation** is a ``mio::TimeSeries`` containing the sizes of each compartment at each time point. A 
simple table can be printed using the ``print_table()`` function of the ``TimeSeries`` class. The output of the 
**FlowSimulation** additionally contains the flows between compartments at each time point. The compartment sizes can 
be printed with ``result[0].print_table()`` and the flows with ``result[1].print_table()``. 
As adaptive step size methods are used by default, the output will not be available on equidistant time points like `dt`
or days. To obtain outputs on days or user-defined time points, you can interpolate the results to days or
any other series of times points with ``mio::interpolate_simulation_result()``.


Visualization
-------------

To visualize the results of a simulation, you can use the Python package :doc:`m-plot <../python/m-plot>`
and its documentation.


List of models
--------------

.. toctree::
    :titlesonly:
    
    models/osir
    models/oseir
    models/oseair
    models/osecir
    models/osecirvvs
    models/osecirts
    models/omseirs4
