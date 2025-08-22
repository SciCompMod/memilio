Linear Chain Trick model
=========================

MEmilio implements a SECIR-type model utilizing the Linear Chain Trick (LCT). This is a generalization of simple ODE-based models and allows for Erlang distributed stay times in the compartments by introducing subcompartments. Note that the resulting system is still described by ODEs. The LCT-SECIR model can be stratified by age groups or other sociodemographic factors.

In the following, we present the general structure of the LCT model. The particular model documentation with examples is linked at the bottom of this page.

An overview of nonstandard but often used data types can be found under :doc:`data_types`.


Infection states
----------------

The model contains a list of **InfectionState**s that define particular features of the subpopulations in the particular state.

.. code-block:: RST

    `State1`
    `State2`
    `...`

To make use of the LCT, we additionally need to define the numbers of subcompartments for each **InfectionState**.

.. code-block:: RST
    `Number of subcompartments of State1`
    `Number of subcompartments of State2`
    `...`

The model is implemented as **CompartmentalModel**.



Sociodemographic stratification
-------------------------------

For this model, the population can be stratified by one sociodemographic dimension. This dimension can be used 
for age groups but for other interpretations. 


Parameters
----------

The parameters of the model are defined as structs and are combined in a class ``ParameterSet<Param1, Param2, ...>``.
We use different types of parameters to represent epidemiological parameters such as the mean stay times in a 
compartment or the contact rates between different age groups. Most model parameters are constants that describe 
pathogen-specific characteristics (possibly resolved by sociodemographic groups) and are represented by a vector with a value for each sociodemographic group. 
To model different contact rates between different sociodemographic groups, we use a parameter denoted **ContactPatterns** of type **UncertainContactMatrix**. 
The **UncertainContactMatrix** contains anmarbitrary large set of contact matrices which can represent the different contact locations in the model like 
schools, workplaces, or homes. The matrices can be loaded or stored in the particular example.
In the **ContactPatterns**, each matrix element stores baseline contact rates :math:`c_{i,j}` between sociodemographic group :math:`i` to group :math:`j`. 
The dimension of the matrix is automatically defined by the model initialization and is reduced to one value if no stratification is used. 
The values can be adjusted during the simulation, e.g., through implementing nonpharmaceutical interventions, 
see the section on :ref:`Nonpharmaceutical Interventions LCT`.
Parameters can be accessed via ``model.parameters.get<Param<double>>()`` and set via either 
``model.parameters.get<Param<double>>() = value`` or ``model.parameters.set<Param<double>>(value)``. 


Initial conditions
------------------

The initial conditions of the model are represented by a class **LctPopulations** that gives the number of individuals in each sociodemographic group and each subcompartment for each **InfectionState**. For more details, see :doc:`Model Creation <lct_creation>`. Before the simulation, the initial conditions for each **InfectionState** and sociodemographic group must be set.


.. _Nonpharmaceutical Interventions LCT:
Nonpharmaceutical interventions
-------------------------------

Contact rates can be adjusted during the simulation to model nonpharmaceutical interventions (NPIs) such as lockdowns, school closures, or social distancing. This is done by adding **Damping**\s to the **ContactPatterns** of the model. A **Damping** is defined by a time point at which the intervention starts and a matrix of the same size as the **ContactMatrix**. 
While in many cases, the reduction matrix is given by a constant matrix with factor :math:`r`, also 
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

Once the model is set up, one can run a simple simulation from time ``t0`` to ``tmax`` with initial step size ``dt`` using the ``mio::simulate()`` function. This will run a simulation of type **Simulation** that saves the sizes of each subcompartment over time. 
You can run a simulation using either fixed or adaptive integration schemes with an absolute or relative tolerance. By default, the simulation uses an adaptive solution scheme of the boost library with an absolute tolerance of 1e-10 and a relative tolerance of 1e-5. For more details on the possible integration schemes, see <numerical integrator stuff>.


Output
------

The output of the **Simulation** is a ``mio::TimeSeries`` containing the sizes of each subcompartment at each time point. 
To obtain a result with respect to the compartments, the subcompartments can be accumulated via the function 
``calculate_compartments()``. A simple table can be printed using the ``print_table()`` function of the 
``mio::TimeSeries`` class. The compartment sizes can be printed with ``model.calculate_compartments(result).print_table()``. 
As adaptive step size methods are used by default, the output will not be available on equidistant time points like `dt` or days. To obtain outputs on days or user-defined time points, you can interpolate the results to days or any other series of times points with ``mio::interpolate_simulation_result()``.


Visualization
-------------

To visualize the results of a simulation, you can use the Python package :doc:`m-plot <../python/m-plot>` and its documentation.


List of models
-----------------------

.. toctree::
    :titlesonly:
    
    models/lsecir
    models/lsecir2d
