IDE models
==========

MEmilio implements two models based on integro-differential equations (IDEs) with different infection states. IDE-based 
models are a generalization of ODE-based models. Whereas ODE-based models assume an exponential distribution regarding 
the time spent in an infection state, IDE-based models allow for arbitrary distributions.

In the following, we present the general structure of the IDE-SECIR model. 
The IDE-SEIR model follows a different structure, an introduction and a detailed example are provided in a separate 
section below. 

Infection States
----------------

Every model contains a list of **InfectionState**\s that define particular features of the subpopulations in the particular state.

.. code-block:: RST

    `State1`
    `State2`
    `...`


Infection State Transitions
---------------------------

Additionally, we define **InfectionTransition**\s that define the possible transitions between the **InfectionState**\s.
When solving the model, we obtain results for the **InfectionTransition**\s as well.

.. code-block:: RST

    `Transition1`
    `Transition2`
    `...`


Sociodemographic Stratification
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
use a parameter denoted **ContactPatterns** of type **UncertainContactMatrix**. The **UncertainContactMatrix** contains
a set of contact matrices of arbitrary length and which can represent the different contact locations in the model like 
schools, workplaces, or homes. The matrices can be loaded or stored in the particular example.

In the **ContactPatterns**, each matrix element stores baseline contact rates :math:`c_{i,j}` between sociodemographic 
group :math:`i` and group :math:`j`. The dimension of the matrix is automatically defined by the model initiation and it is reduced 
to one value if no stratification is used. The values can be adjusted during the simulation, e.g., through implementing 
nonpharmaceutical interventions, see the section on :ref:`Nonpharmaceutical Interventions`. 

An important feature of our IDE-based model is that we can choose the transition distributions in a flexible way. The 
default distribution is a smoother cosine function as it provides good testing qualities. For more realistic simulations, 
MEmilio provides the possibility to use exponential, gamma or lognormal distributions within the model.
Practically, one first needs to create an object of a class that is derived from the class ``StateAgeFunction``, 
e.g. ``SmootherCosine``. Any class that is derived from ``StateAgeFunction`` can be inserted into a 
``StateAgeFunctionWrapper`` object that is then passed to the model.

Parameters can get accessed via ``model.parameters.get<Param<ScalarType>>()`` and set via either 
``model.parameters.get<Param<ScalarType>>() = value``.



Initial Conditions
------------------

The initial conditions consist of a **TimeSeries** that contains the flows per time step for some time interval before 
the simulation start. The length of this time interval depends on the chosen transition distributions and takes into 
account how long there is a relevant fraction of individuals remaining in the compartments. For more information, see 
the documentation of **StateAgeFunction**. Note that the last time point of the initial **TimeSeries** determines the 
start time of the simulation. 


.. _Nonpharmaceutical Interventions:
Nonpharmaceutical Interventions
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
-----------------------

.. toctree::
    :titlesonly:
    
    models/iseir
    models/isecir



IDE-SEIR model
---------------

Introduction
~~~~~~~~~~~~~
The four compartments 

- `Susceptible` (:math:`S`), may become exposed at any time
- `Exposed` (:math:`E`), becomes infected after some time
- `Infected` (:math:`I`), will recover after some time
- `Recovered` (:math:`R`)

are used to simulate the spread of the disease. 

Simulation
~~~~~~~~~~~

The simulation runs in discrete time steps using a trapezoidal rule. The model and the numerical scheme is based on the paper `"Modeling infectious diseases using integro-differential equations: Optimal
control strategies for policy decisions and Applications in COVID-19" by Keimer and Pflug, 2020 <http://dx.doi.org/10.13140/RG.2.2.10845.44000>`_. 

For a detailed description and application of the model, see:

Plötzke L (2021) Modellierung epidemischer Infektionskrankheiten auf der Basis von gewöhnlichen und Integro-Differentialgleichungen. Bachelor thesis, University of Cologne. https://elib.dlr.de/143504/

How to: Set up and run a simulation of the IDE-SEIR model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


Here, we set the contact matrix used in the simulation. One can define multiple matrices for different locations. The size of each of these matrices is defined by the number of age groups. 
Below, we use only one contact matrix for one location. As we only consider one age group in our example, we set the corresponding contact rate to :math:`10`.

.. code-block:: cpp

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));

To simulate the implementation of nonpharmaceutical interventions, we add dampings to the contact rate. Here, we apply a damping of :math:`0.7` after :math:`10` days, meaning that the contact rate is reduced to :math:`30%` of the initial value.  

.. code-block:: cpp

    contact_matrix[0].add_damping(0.7, mio::SimulationTime(10.));
    model.parameters.get<mio::iseir::ContactFrequency<double>>() = mio::UncertainContactMatrix<double>(contact_matrix);

After defining :math:`t_{\max}`, we can simulate, which means that we calculate the value for the compartment :math:`S`.

.. code-block:: cpp

    int tmax  = 15;
    model.simulate(tmax);

The values of the remaining compartments :math:`E`, :math:`I` and :math:`R` are calculated using the parameters ``LatencyTime`` and ``InfectiousTime`` and obtain a time series containing the values of all compartments. 

.. code-block:: cpp

    auto result = model.calculate_EIR();

Finally, we can print our results. 

.. code-block:: cpp

    result.print_table({"S", "E", "I", "R"});


