Tutorial: Usage of Python bindings
==================================

This tutorial should give an overview of how to use the
currently available functions and models of the Python bindings.
For expanding the bindings with new models look into the section 
:doc:`Model Creation <model_creation>`.

Generally, the package is following the structure of the main C++
library to make it easy to understand and comparable, while introducing
changes to create a more pythonic interface. We will follow the example of an
ODE SEIR model.

Define infectious disease model
--------------------------------

Following is a comparison of the model initialization between
Python and C++ to better understand the differences of both interfaces.

.. grid:: 1 1 2 2

   .. grid-item::
      
      Python:

      .. code-block:: python
         
         import memilio.simulation.oseir as oseir

         num_groups = 1
         model = oseir.Model(num_groups)

   .. grid-item::
      
      C++:

      .. code-block:: cpp


         #include "ode_seir/model.h"

         int num_groups = 1;
         mio::oseir::Model<double> model(num_groups);

Initialize parameters
---------------------

Next, the parameters should be defined of the scenario you want to model. The model has the 
possibility to incorporate age stratification, which leads to many of the parameters having 
values for each age group. The AgeGroup is used for indexing. For the example model
with a single age group defined by num_groups the parameter definitions follow as

.. code-block:: python

   import memilio.simulation as mio

   A0 = mio.AgeGroup(0)

   # Compartment transition duration
   model.parameters.TimeExposed[A0] = 5.2
   model.parameters.TimeInfected[A0] = 6.

   # Compartment transition propabilities
   model.parameters.TransmissionProbabilityOnContact[A0] = 1.

``AgeGroup(0)`` defines the first age group. For a model with more than one age group,
we could index the other groups with ``AgeGroup(1)``, ``AgeGroup(2)``, ....

Initial conditions
-------------------

We also need to define the inital states of the population. They are not only divided through an age group,
but also an infection state, such that an additional index of the enum ``InfectionState`` has to be provided.

.. code-block:: python

   total_population = 83_000

   model.populations[A0, oseir.InfectionState.Exposed] = 100
   model.populations[A0, oseir.InfectionState.Infected] = 50
   model.populations[A0, oseir.InfectionState.Recovered] = 10
   model.populations.set_difference_from_total(
      (A0, oseir.InfectionState.Susceptible), total_population)

The function ``model.populations.set_difference_from_total()`` helps by setting the last compartment with
a total population size of the scenario and will take the difference to the combined other compartments
of the age group.

Nonpharmaceutical interventions
-------------------------------

One important topic of interest in infectious disease are nonpharmaceutical interventions, which aim to reduce the spread of the disease.
The model supports to define base contact patters for the mean daily contacts between different age groups. When using multiple contact patterns they get 
accumulated to get overall daily contacts. 

.. code-block:: python

   model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
        (num_groups, num_groups))
   model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.zeros(
        (num_groups, num_groups))

Then, dampings can be added to (partially) reduce the contacts defined by a ``ContactMatrix`` beginning at a time step ``t``. 

.. code-block:: python

   model.parameters.ContactPatterns.cont_freq_mat.add_damping(
        Damping(coeffs=np.r_[0.9], t=30.0, level=0, type=0))

Simulation
----------

Now, the infectious diesease dynamic can be simulated by calling ``simulate()``:

.. code-block:: python

   result = oseir.simulate(t0=0, tmax=60, dt=1, model)

Similar to the MEmilio C++ library, the Python interface provides the option of adjusting the solver.
Currently available are:

* Euler
* RungeKuttaCashKarp54 (default)
* RungeKutta-Fehlberg45

The integrator can be changed as the last parameter of the simulate function.

.. code-block:: python

   integrator = mio.RKIntegratorCore(dt_max=1)
   result = oseir.simulate(0, tmax=60, dt=1, model, integrator)

Output and visualization
-------------------------

The result returned from the simulation is a TimeSeries object containing the number of people per age group in each infection state at each time step.
The TimeSeries provides alot of interfaces to interact with it, but can also be transformed into a multidimensional numpy matrix for a more
pythonic interface.

.. code-block:: python
   
   result_array = result.as_ndarray()

Now you can use the usual data handling options and make us of the easy visualization tools that are part of Python.
Some plotting functions specific to MEmilio and created as part of the project are combined in the :doc:`MEmilio Plot Package <memilio_plot>`.

Additional ressources
---------------------

Further examples are provided at `examples/simulation <https://github.com/SciCompMod/memilio/blob/main/pycode/examples/simulation/>`_. 
They include the usage of a FlowModel, introducing a graph model for regional differences or parameter studies for simulating under uncertainty.



