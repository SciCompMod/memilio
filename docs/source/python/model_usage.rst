Tutorial: Usage of python bindings
==================================

This tutorial should give an overview of how to use the
currently available functions and models of the python bindings.
For expanding the bindings with new models look into the section 
:doc:`Model Creation <model_creation>`.

Generally, the package is following the structure of the main C++
library to make it easy to understand and comparable, while introducing
changes to create a more pythonic interface. We will follow the example of an
ODE SEIR model starting with a comparison of the model initialization between
Python and C++.

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

         mio::oseir::Model<ScalarType> model(1);


Next, the parameters should be defined of the scenario you want to model. The model has the 
possibility to incorporate age stratification, which leads to many of the parameters to have 
values for each age group. The AgeGroup is used for indexing, which leads for the example model
with a single group defined by num_groups to the parameter definitions.

.. code-block:: python

   A0 = AgeGroup(0)

   # Compartment transition duration
   model.parameters.TimeExposed[A0] = 5.2
   model.parameters.TimeInfected[A0] = 6.

   # Compartment transition propabilities
   model.parameters.TransmissionProbabilityOnContact[A0] = 1.

AgeGroup(0) defines the first age group. For a model with more than one age group,
we could index the other groups with AgeGroup(1), AgeGroup(2), ....

We also need to define the inital states of the population. As they are not only divided through an age group,
but also an infection state, we need to add an index of the enum InfectionState.

.. code-block:: python

   total_population = 83_000

   model.populations[A0, InfectionState.Exposed] = 100
   model.populations[A0, InfectionState.Infected] = 50
   model.populations[A0, InfectionState.Recovered] = 10
   model.populations.set_difference_from_total(
      (A0, InfectionState.Susceptible), total_population)

The function model.populations.set_difference_from_total hepls by setting the last compartment with
a total population size of the scenario and will take the difference to the combined other compartments
of the age group.


Now we could simulate with:

.. code-block:: python

   result = simulate(0, days, dt, model)

Similar to the MEmilio C++ library, the python interface provides the option of adjusting the solver.
Currently available:
-
-

.. code-block:: python

   result = simulate(0, days, dt, model)

parameter setting with spatial resolution

simulation
redefining the solver

Working with the output

Expanding to graph model

More examples:
- `examples/ode_secir.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secir.cpp>`_
- `examples/ode_secir_ageres.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secir_ageres.cpp>`_
- `examples/ode_secir_parameter_study.cpp <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/ode_secir_parameter_study.cpp>`_


Lastly Limitations with introduction of model creation

