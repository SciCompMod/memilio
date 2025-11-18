How to: Usage of Python bindings
==================================

This tutorial should give an overview of how to use the
currently available functions and models of the Python bindings.
For expanding the bindings with new models look into the section 
:doc:`Model Creation <m-simulation_expanding_bindings>`.

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

One important topic of interest in infectious disease dynamics are nonpharmaceutical interventions (NPIs), which aim to mitigate the dynamics of the disease spread. 
In our models, NPIs are implemented through dampings to the contact matrix. These dampings reduce the contact rates between different groups to simulate interventions.

The basic contact matrix is defined through its baseline `ContactPatterns` (and a potential minimum pattern which does not need to be set and which is not used by default).
For a model with a single age group, the contact matrix is simply value, a 1x1 matrix and `num_groups=1`. For `num_groups>1`, the contact matrix is a square matrix of size `num_groups x num_groups`,
where each entry represents the contact frequency between two age groups.

.. code-block:: python

   model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
        (num_groups, num_groups))
   model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.zeros(
        (num_groups, num_groups))

Then, dampings can be added to (partially) reduce the contacts defined by a ``ContactMatrix`` beginning at a time step ``t=30``. 

.. code-block:: python

   model.parameters.ContactPatterns.cont_freq_mat.add_damping(
        Damping(coeffs=np.r_[0.9], t=30.0, level=0, type=0))


If a minimum pattern is set, the contact reduction will relatively reduced the difference between the baseline and minimum patterns.

Several models also supports dynamic NPIs based on epidemic thresholds. These are implemented in the model specific **Simulation** class and are automatically triggered based on predefined criteria, such as the percentage of infected individuals in the population.

For more complex scenarios, such as real-world lockdown modeling, you can implement detailed NPIs with location-specific dampings (e.g., home, school, work, other). 

Example for defining different contact locations:

.. code-block:: python

   // Define different contact locations
   class Location(Enum):
      """ """
      Home = 0
      School = 1
      Work = 2
      Other = 3

    
   // Load contact location matrices for Germany
   contact_matrices = mio.ContactMatrixGroup(
      len(list(Location)), self.num_groups)
   locations = ["home", "school_pf_eig", "work", "other"]

   for i, location in enumerate(locations):
      baseline_file = os.path.join(
            self.data_dir, "Germany", "contacts", "baseline_" + location + ".txt")
      contact_matrices[i] = mio.ContactMatrix(
            mio.read_mobility_plain(baseline_file),
            np.zeros((self.num_groups, self.num_groups))
      )
   model.parameters.ContactPatterns.cont_freq_mat = contact_matrices

You can create intervention types that target specific locations with different intensities:

.. code-block:: python

    // Different types of NPI
   class Intervention(Enum):
      """ """
      Home = 0
      SchoolClosure = 1
      HomeOffice = 2
      GatheringBanFacilitiesClosure = 3
      PhysicalDistanceAndMasks = 4
      SeniorAwareness = 5
    
    // Different levels of NPI
   class InterventionLevel(Enum):
      """ """
      Main = 0
      PhysicalDistanceAndMasks = 1
      SeniorAwareness = 2
      Holidays = 3

A complex lockdown scenario with multiple interventions starting on a specific date can be implemented via:

.. code-block:: python

   start_spring_date = datetime.date(2020, 3, 18)

   if start_spring_date < end_date:
      start_spring = (start_spring_date - self.start_date).days
      dampings.append(contacts_at_home(start_spring, 0.6, 0.8))
      dampings.append(school_closure(start_spring, 1.0, 1.0))
      dampings.append(home_office(start_spring, 0.2, 0.3))
      dampings.append(social_events(start_spring, 0.6, 0.8))
      dampings.append(social_events_work(start_spring, 0.1, 0.2))
      dampings.append(physical_distancing_home_school(start_spring, 0.4, 0.6))
      dampings.append(physical_distancing_work_other(start_spring, 0.4, 0.6))

A more advances structure to automatically activate interventions based on threshold criteria is given by **DynamicNPIs**.
Dynamic NPIs can be configured to trigger when the number of symptomatic infected individuals exceeds a certain relative threshold in the population. 
In contrast to static NPIs which are active as long as no other NPI gets implemented, dynamic NPIs are checked at regular intervals and get 
activated for a defined duration when the threshold is exceeded. As above, different dampings `dampings` can be assigned to different contact locations
and are then triggered all at once the threshold is exceeded.
The following example shows how to set up dynamic NPIs based on the number of 200 symptomatic infected individuals per 100,000 population. 
It will be active for at least 14 days and checked every 3 days. If the last check after day 14 is negative, the NPI will be deactivated.

.. code-block:: python

    // Configure dynamic NPIs with thresholds
   dynamic_npis = params.DynamicNPIsInfectedSymptoms
   dampings = []
   # increased from [0.4, 0.6] in Nov
   dampings.append(contacts_at_home(0, 0.6, 0.8))
   dampings.append(school_closure(0, 0.25, 0.25))  # see paper
   dampings.append(home_office(0, 0.2, 0.3))
   dampings.append(social_events(0, 0.6, 0.8))
   dampings.append(social_events_work(0, 0.1, 0.2))
   dampings.append(physical_distancing_home_school(0, 0.6, 0.8))
   dampings.append(physical_distancing_work_other(0, 0.6, 0.8))
   dampings.append(senior_awareness(0, 0.0, 0.0))

   dynamic_npis.interval = 3.0
   dynamic_npis.duration = 14.0
   dynamic_npis.base_value = 100000
   dynamic_npis.set_threshold(200.0, dampings)

            

Simulation
----------

Now, the infectious diesease dynamics can be simulated by calling ``simulate()``:

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

Now you can use the usual data handling options and make use of the easy visualization tools that are part of Python.
Some plotting functions specific to MEmilio and created as part of the project are combined in the :doc:`MEmilio Plot Package <m-plot>`.

Additional resources
---------------------

Further examples are provided at `examples/simulation <https://github.com/SciCompMod/memilio/blob/main/pycode/examples/simulation/>`_. 
They include the usage of a FlowModel, introducing a graph model for regional differences or parameter studies for simulating under uncertainty.



