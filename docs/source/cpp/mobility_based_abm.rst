Agent-based model
=================

This module models and simulates the epidemic using an agent-based model (*ABM*) approach. Unlike the compartmental models that use a system of ODEs, this model simulates
the spread of COVID-19 in a population with discrete persons (the agents) moving throughout locations in the
model and interacting with (infecting) each other. For a detailed overview of the ABM, see:

- Kerkmann D, Korf S, et al. Agent-based modeling for realistic reproduction of human mobility and contact behavior to evaluate test and isolation strategies in epidemic infectious disease spread. https://doi.org/10.48550/arXiv.2410.08050

Introduction
-----------

The model is implemented in multiple classes and headers located in the ``/cpp/models/abm/`` directory. The core classes and their locations are:

- ``Person`` (person.h): Represents individual agents in the simulation
- ``Infection`` (infection.h): Manages infection dynamics and disease progression 
- ``Location`` (location.h): Defines places where agents interact
- ``Model`` (model.h): Coordinates all components of the simulation
- ``Simulation`` (simulation.h): Executes the simulation logic

The following sections outline the major features of the agent-based model.

Structure
~~~~~~~~~

The model consists of the following major classes:

1. **Person**: Represents an agent of the model. A person has an ID and a list with their assigned locations (i.e. the locations they visit during the simulation). They can perform
   tests and wear masks. Every person has lists with past and current infections and vaccinations.
   
2. **Infection**: Collection of all information about a person's infection, i.e. infectiousness, infection course,
   virus variant. The infection course is drawn stochastically from the infection states that are similar to the
   compartments of the SECIR model.
   
3. **Location**: Represents places in the model where people meet and interact, e.g. home, school, work, social event
   sites. A location can be split into cells to model parts of a location, like classrooms in a school. Some infection
   parameters are location-specific and one can activate NPIs like mandatory masks or tests to enter the location.
   
4. **Model**: Collection of all persons and locations. It also holds information about the testing strategy of the
   simulation and holds the rules for the mobility phase.
   
5. **Simulation**: Runs the simulation and stores results.

Disease progression
~~~~~~~~~~~~~~~~~~

The ABM implements a detailed disease progression model that captures the full course of an infection from exposure to resolution. The disease progression is modeled through the ``Infection`` class, which contains:

1. **Infection States**: Similar to the SECIR model, an infected person progresses through states defined in ``infection_state.h``:

   * **Susceptible**: Initial state before infection
   * **Exposed**: Infected but not yet infectious
   * **InfectedNoSymptoms**: Infectious but without symptoms
   * **InfectedSymptoms**: Showing symptoms but not severe
   * **InfectedSevere**: Severe infection requiring hospitalization
   * **InfectedCritical**: Critical infection requiring ICU
   * **Recovered**: Recovered from infection with immunity
   * **Dead**: Deceased due to infection

2. **Viral Load Dynamics**: The model implements realistic viral load curves based on scientific data:

   * **Incline Phase**: Rapid increase in viral concentration
   * **Peak**: Maximum viral load
   * **Decline Phase**: Gradual decrease until clearance
   
3. **Infectiousness**: The probability of transmitting the virus depends on viral load through an invlogit function.

4. **Stochastic Transitions**: Progression between states is stochastic, with age-dependent probabilities:

   * The duration in each state is drawn from distributions in the model parameters
   * Prior immunity (from vaccination or previous infection) affects:

     * Viral load (reduced peak)
     * Severity progression (reduced probability of severe outcomes)
     * Duration of infectious period
   
5. **Infection Course**: The infection course is determined by:

   * Age group of the person
   * Virus variant
   * Protection status (prior immunity)
   * Random factors (individual variation)

Data collection
~~~~~~~~~~~~~~~~~~

The ABM simulation can collect data through the ``History`` object, which allows for flexible data logging. This is particularly 
useful for analyzing results after the simulation has completed. There are multiple types of data that can be collected:

1. **Time Series Data**: Track how infection states change over time
   
2. **Location-specific Data**: Monitor occupancy or infection rates at specific locations

3. **Person-specific Data**: Follow individual movement patterns or infection trajectories

The examples demonstrate two approaches:

.. code-block:: cpp

   // Basic time series tracking of infection states
   mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
       Eigen::Index(mio::abm::InfectionState::Count)};
   
   // More complex logging with multiple data types
   mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds> history;
   
   // Run simulation with history object
   sim.advance(tmax, history);

Interventions
~~~~~~~~~~~~~~~~~~

The ABM supports various interventions that can be applied at specific time points, such as:

1. **Capacity Restrictions**: Limit the number of people at locations

2. **Testing Regimes and Quarantines**: Implement regular testing at specific locations and resulting quarantines at home

3. **Lockdowns**: Restrict movement between locations

Simulation
----------

The simulation runs in discrete time steps. Each step has two phases, an **interaction phase** and a **mobility phase**.
After these two phases the disease can progress and the simulation time is increased by one step.

Interaction phase
~~~~~~~~~~~~~~~~~~~

In this phase, each person interacts with the other persons at the same location. This interaction determines the
transmission of the disease. A susceptible person can become infected by contact with an infected person. The probability
of infection depends on a multitude of factors, such as the viral load and infectiousness of the infected and the immunity
level of the susceptible person.

Mobility phase
~~~~~~~~~~~~~~~~~~

During the mobility phase, each person may change their location. Mobility follow
`rules <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/mobility_rules.cpp>`_, considering the current location, time of day, and properties of the person (e.g. age).
Some location changes are deterministic and regular (e.g. going to work), while others are random (e.g. going shopping or to a
social event in the evening/on the weekend). When agents are infected, they are quarantined and cannot change their location.
You can restrict some mobility rules by allowing only a proportion of people to enter specific locations.

Another way of mobility is using trips. A trip consists of the ID of the person that performs this trip, a time point when this trip is performed, and the destination.
At the beginning of the simulation, a list with all trips is initialized and followed during the simulation. There can be different
trips on the weekend than during the week, but other than that, the agents do the same trips every day. As before, agents that are
in quarantine or in the hospital cannot change their location.

How to
-----------

This section gives an introduction to how to use the ABM and set up your own simulation. For a quick overview, you can find a full
example in the `ABM minimal example <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/abm_minimal.cpp>`_ and a more detailed Doxygen documentation
`here <https://scicompmod.github.io/memilio/documentation/index.html>`_. For a guide on installation and running the simulations and
examples, see this `README <https://github.com/SciCompMod/memilio/blob/main/cpp/README.md>`_.

Every person in the ABM belongs to an AgeGroup, which we can define as follows:

.. code-block:: cpp

   size_t num_age_groups         = 4;
   const auto age_group_0_to_4   = mio::AgeGroup(0);
   const auto age_group_5_to_14  = mio::AgeGroup(1);
   ...                           = ...

Note that every age group has to have values strictly smaller than the number of age groups ``num_age_groups``.
With this number we create an empty model:

.. code-block:: cpp

   auto model = mio::abm::Model(num_age_groups);

We can set several general parameters, which you can find `here <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/parameters.h>`_. Here is an example where we set the
duration of the incubation period to 4 days:

.. code-block:: cpp

   model.parameters.get<mio::abm::IncubationPeriod>() = 4.;

Locations and persons
~~~~~~~~~~~~~~~~~~~~~

To add a location to the model, we have to specify the kind of location:

.. code-block:: cpp

   auto home = model.add_location(mio::abm::LocationType::Home);

People are added with an age. Then we have to assign them, so the model knows they can travel to this location:

.. code-block:: cpp

   auto person = model.add_person(home, age_group_0_to_4);
   person.set_assigned_location(home);

For more complex location configurations, the model allows setting location-specific parameters:

.. code-block:: cpp

   // Add one social event with 5 maximum contacts
   auto event = model.add_location(mio::abm::LocationType::SocialEvent);
   model.get_location(event).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
   
   // Increase aerosol transmission for all locations
   model.parameters.get<mio::abm::AerosolTransmissionRates>() = 10.0;
   
   // Increase contact rate for specific age groups at work
   model.get_location(work)
       .get_infection_parameters()
       .get<mio::abm::ContactRates>()[{age_group_15_to_34, age_group_15_to_34}] = 10.0;

Households
~~~~~~~~~~

For adding more people to the model, we can create households. A Household holds a vector of HouseholdMembers, which in turn
hold a weighted distribution, such that we can randomly draw the age of each Person belonging to the Household. To manage
multiple Households of the same type, we can use a HouseholdGroup.
In our example, we categorize individuals into two groups: children and parents.

.. code-block:: cpp

   auto child = mio::abm::HouseholdMember(num_age_groups);
   child.set_age_weight(age_group_0_to_4, 1);
   child.set_age_weight(age_group_5_to_14, 1);

   auto parent = mio::abm::HouseholdMember(num_age_groups);
   parent.set_age_weight(age_group_15_to_34, 1);
   parent.set_age_weight(age_group_35_to_59, 1);

   // Two-person household with one parent and one child.
   auto twoPersonHousehold_group = mio::abm::HouseholdGroup();
   auto twoPersonHousehold_full  = mio::abm::Household();
   twoPersonHousehold_full.add_members(child, 1);
   twoPersonHousehold_full.add_members(parent, 1);
   twoPersonHousehold_group.add_households(twoPersonHousehold_full, n_households);
   add_household_group_to_model(model, twoPersonHousehold_group);

Testing strategies
~~~~~~~~~~~~~~~~~

During the simulation, people can get tested, and we have to specify the scheme for that:

.. code-block:: cpp

   auto validity_period       = mio::abm::days(1);
   auto probability           = 0.5;
   auto start_date            = mio::abm::TimePoint(0);
   auto end_date              = mio::abm::TimePoint(0) + mio::abm::days(30);
   auto test_type             = mio::abm::TestType::Antigen;
   auto test_parameters       = model.parameters.get<mio::abm::TestData>()[test_type];
   auto testing_criteria_work = mio::abm::TestingCriteria();
   auto testing_scheme_work   = mio::abm::TestingScheme(testing_criteria_work, validity_period, 
                                                     start_date, end_date,
                                                     test_parameters, probability);
   model.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme_work);

Initializing infections
~~~~~~~~~~~~~~~~~~~~~~

For some infections to happen during the simulation, we have to initialize people with infections:

.. code-block:: cpp

   // Assign infection state to each person randomly with specific distribution
   std::vector<double> infection_distribution{0.5, 0.3, 0.05, 0.05, 0.05, 0.05, 0.0, 0.0};
   for (auto& person : model.get_persons()) {
       mio::abm::InfectionState infection_state = mio::abm::InfectionState(
           mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), infection_distribution));
       auto rng = mio::abm::PersonalRandomNumberGenerator(person);
       if (infection_state != mio::abm::InfectionState::Susceptible) {
           person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, 
                                                       person.get_age(),
                                                       model.parameters, start_date, infection_state));
       }
   }

Running the simulation
~~~~~~~~~~~~~~~~~~~~~

Finally, we run the simulation:

.. code-block:: cpp

   auto t0   = mio::abm::TimePoint(0);
   auto tmax = t0 + mio::abm::days(30);
   auto sim  = mio::abm::Simulation(t0, std::move(model));
   
   // Simple simulation without data collection
   sim.advance(tmax);

Alternatively, if we want to track things in the simulation, we need to set up a
`history <https://github.com/SciCompMod/memilio/blob/main/cpp/memilio/io/README.md#the-history-object>`_, for example, to track all the Infection states of each simulation step.

.. code-block:: cpp

   mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> history{
       Eigen::Index(mio::abm::InfectionState::Count)};

Then we can run the simulation with the history object and access the data through ``get_log()``:

.. code-block:: cpp

   sim.advance(tmax, history);
   auto log = history.get_log();

Finally, we can print the data to a text file:

.. code-block:: cpp

   std::ofstream outfile("abm_minimal.txt");
   std::get<0>(log).print_table({"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4, outfile);
   std::cout << "Results written to abm_minimal.txt" << std::endl;

