.. include:: ../literature.rst

Agent-based model
=================

This module models and simulates the epidemic using an agent-based model (*ABM*) approach. Unlike the compartmental models that use a system of ODEs, this model simulates
the spread of an epidemic in a population with discrete persons (the agents) moving throughout locations in the
model and interacting with (infecting) each other. The model aims to capture the dynamics of disease spread in a realistic way, considering individual behaviors, mobility patterns, and disease progression.
Thus, it allows for a more detailed and flexible representation of the epidemic dynamics, including the effects of interventions and individual-level variations.
Consequently, the ABM is computationally more expensive than the ODE models, but it provides a more granular view of the epidemic spread.
At the same time, the ABM has to be parametrized carefully to ensure that the simulation results are realistic and meaningful.

In the following, we present more details of the agent-based model, including code examples. 
An overview of nonstandard but often used data types can be found under :doc:`data_types`.

Structure
~~~~~~~~~

The model is implemented in multiple classes. Source and header files are located in the `/cpp/models/abm/ <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/>`_ directory. While many files contain supporting implementation or additional features, it is important to understand the main workflow and core classes.
The core classes and their locations are:

- ``Simulation`` (`simulation.h <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/simulation.h>`_): Runs the simulation and stores results. The model is evolved in discrete time steps of the same size which can be chosen by the user.
- ``Model`` (`model.h <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/model.h>`_): Collection of all persons, locations and used parameters. Is initialized with the number of age groups that are considered. It also holds information about the testing strategy of the simulation and holds the rules for the mobility phase.
- ``Model Functions`` (`model_functions.h <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/model_functions.h>`_): A collection of model functions that mostly cover interaction of agents at locations. These functions are called in the model evolve function.
- ``Parameters`` (`parameters.h <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/parameters.h>`_): Collection of all parameters used in the model.
- ``Location`` (`location.h <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/location.h>`_): Represents places in the model where people meet and interact, e.g. home, school, work, social event sites. Along their type, locations contain an ID and a geographical location (longitude and latitude). A location can be split into cells to model parts of a location, like classrooms in a school. Some infection parameters are location-specific and can be set per location. Mandatory masks can be activated to simulate a mask obligation intervention.
- ``Person`` (`person.h <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/person.h>`_): Represents an agent of the model. A person has an ID, is associated with an age group and has a list with their assigned locations (i.e. the locations they can visit during the simulation). Every person has lists with past and current infections as well as vaccinations. Further, more information on the personal behavior and test results is available.
- ``Infection`` (`infection.h <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/infection.h>`_): Collection of all information about a person’s infection, i.e. infectiousness, infection course and symptoms, virus variant. The infection course is drawn stochastically from the infection states that are similar to the compartments of the SECIR model and is explained in detail below.


Disease progression
~~~~~~~~~~~~~~~~~~~

The ABM implements a detailed disease progression model that captures the full course of an infection from exposure to resolution. The disease progression is modeled through the ``Infection`` class, which contains:

1. **Infection States**: Similar to the aggregated models (see, e.g., :doc:`equation based models<ode>`), an infected person progresses through states defined in `infection_state.h <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/infection_state.h>`_:

   * **Susceptible**: Initial state before infection
   * **Exposed**: Infected but not yet infectious
   * **InfectedNoSymptoms**: Infectious but without symptoms
   * **InfectedSymptoms**: Showing symptoms but not severe
   * **InfectedSevere**: Severe infection requiring hospitalization
   * **InfectedCritical**: Critical infection requiring ICU
   * **Recovered**: Recovered from infection with immunity
   * **Dead**: Deceased due to infection

Stochastic Transitions:
   Agents traverse the infection states from Susceptible to Recovered or Dead. Recovery is possible from every Infected State (NoSymptoms, Symptoms, Severe, Critical) and Dead is possible from InfectedSevere and InfectedCritical.
   Right now, agents that have recovered from an infection cannot be infected again. This is supposed to change in the future, allowing for reinfections.
   Progression between states is stochastic, with age-dependent probabilities. The duration in each state is drawn from distributions.
   Note that vaccination is not implemented through infection states, but separately and has an impact on the disease transmission and progression (compare with 4. **Dependencies** below).

2. **Viral Load Dynamics**: The model implements realistic viral load curves based on scientific data. The logarithmic viral load is modeled as a function of time since infection, with three phases:

   * **Incline Phase**: Rapid increase in viral concentration
   * **Peak**: Maximum viral load
   * **Decline Phase**: Gradual decrease until clearance
   
   From the individual viral load, the viral shed is calculated, which represents the amount of virus that is shed by an infected person and can potentially infect others.
   This quantity is used in the interaction phase to determine the probability of transmission during interactions.

3. **Dependencies**: Both infection state and viral load parameters depend on:

   * Age group of the person
   * Virus variant
   * Protection status (prior immunity)
   * Random factors (individual variation)

   In particular, prior immunity (from vaccination or previous infection) affects:

     * Viral load (reduced peak)
     * Severity progression (reduced probability of severe outcomes)
     * Duration of infectious period
     * Probability of being infected (again)

4. **Disease spread**: During interactions, agents can infect each other. The viral shed is used in combination with further personal information and contact details to feed into a stochastic process that determines if the virus is transmitted and a new agent becomes infected. The chosen time step of the model has no impact on the expected amount of transmissions per time.

For details on the mathematical modeling of viral shed and disease spread, we refer to 

- |Agent-based_modeling_for|

Data extraction
~~~~~~~~~~~~~~~
The ABM simulation can collect and extract data through the ``History`` object, which allows for flexible data logging and writing.
A collection of often used loggers and writers is available in `common_abm_loggers.h <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/common_abm_loggers.h>`_, but users can define their own loggers and writers to satisfy their individual needs.
This is particularly useful for analyzing results after the simulation has completed. There are multiple types of data that can be collected:

1. **Time Series Data**: Track how infection states change over time
   
2. **Location-specific Data**: Monitor occupancy or infection rates at specific locations

3. **Person-specific Data**: Follow individual movement patterns or infection trajectories

However, the user can also define their own loggers and writers to collect custom data. For further information and examples, see :doc:`io`.

Interventions
~~~~~~~~~~~~~

The ABM supports various interventions that can be applied at specific time points, such as:

1. **Capacity Restrictions**: Limit the number of people at locations

2. **Testing Regimes and Quarantines**: Implement regular testing at specific locations and resulting quarantines at home

3. **Lockdowns**: Restrict movement between locations

Examples for usage can be found below.

Simulation
----------

The simulation runs in discrete time steps. Each step has two phases, an **interaction phase** and a **mobility phase**.
After these two phases, the disease can progress and the simulation time is increased by one step.

Interaction phase
~~~~~~~~~~~~~~~~~

In this phase, each person interacts with the other persons at the same location. This interaction determines the
transmission of the disease. A susceptible person can become infected by contact with an infected person. The probability
of infection depends on a multitude of factors, such as the viral load and infectiousness of the infected and the immunity
level of the susceptible person.

Mobility phase
~~~~~~~~~~~~~~

During the mobility phase, each person may change their location.

The available location types defined in `location_type.h <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/location_type.h>`_ are:

   * **Home**: Home location of a person
   * **School**: School location for children
   * **Work**: Workplace for adults
   * **SocialEvent**: Locations for social gatherings (e.g., parties, events)
   * **BasicsShop**: Basic shop for essential goods (e.g., grocery store)
   * **Hospital**: Hospital for severely infected persons
   * **ICU**: Intensive Care Unit for critical patients
   * **Cemetery**: Exists once per model and is used as a final resting place for deceased persons

A few more types are available, but these are currently not used in the model.

The model supports two ways of mobility:
`Mobility rules <https://github.com/SciCompMod/memilio/blob/main/cpp/models/abm/mobility_rules.cpp>`_, considering the current location, time of day, and properties of the person (e.g. age).
The mobility rules use the assigned locations of the persons. Some location changes are deterministic and regular (e.g. going to work), while others are random (e.g. going shopping or to a
social event in the evening/on the weekend). When agents are infected, they are quarantined and cannot change their location.
You can restrict some mobility rules by allowing only a proportion of people to enter specific locations. We divide the mobility rules into two categories:

1. **Infection-based mobility**: This mobility is based on the infection state of the person. For example, a person in quarantine cannot change their location, and severely or critically infected persons go to the hospital or ICU.
   This mobility is used to model the behavior of people during an epidemic. It consists of the following rules:

   * Going home when quarantined
   * Going to the hospital when severely infected
   * Going to the ICU when critically infected
   * Going to the cemetery when deceased
   * Returning home when recovered

   More severe cases of infection take precedence over less severe cases, meaning for example that a critically infected person goes to the ICU, and does not stay in quarantine at home.

2. **Optional mobility**: This mobility is not based on the infection state of the person. For example, a person can go to a social event or a shop.
While the first category is mandatory, the second category is optional and can be restricted by the user. This allows for modeling different scenarios, such as lockdowns or social distancing measures, or the exclusive usage of trips.
The optional mobility rules consist of:

   * Going to work at work hours
   * Going to school at school hours
   * Going to a social event in the evening or on weekends
   * Going to a shop randomly during the day (except Sunday)

Another way of mobility is using trips. A trip consists of the ID of the person that performs this trip, a time point when this trip is performed, and the destination.
At the beginning of the simulation, a list with all trips is initialized and followed during the simulation. The agents do the same trips every day. As before, agents that are
in quarantine or in the hospital cannot change their location. Trips can be used even for locations that are not the assigned locations for the respective person.


How to
------

This section gives an introduction to how to use the ABM and set up your own simulation. For a quick overview, you can find a full
example in the `ABM minimal example <https://github.com/SciCompMod/memilio/blob/main/cpp/examples/abm_minimal.cpp>`_. For a guide on installation and running the simulations and
examples, see :doc:`installation`.

Every person in the ABM belongs to an AgeGroup, which we can define as follows:

.. code-block:: cpp

   size_t num_age_groups         = 4;
   const auto age_group_0_to_4   = mio::AgeGroup(0);
   const auto age_group_5_to_14  = mio::AgeGroup(1);
   ...                           = ...

Note that every age group has to have arguments strictly smaller than the number of age groups ``num_age_groups``.
With this number we create an empty model:

.. code-block:: cpp

   auto model = mio::abm::Model(num_age_groups);

The model parameters can be set for the whole model or for specific locations. For example, we can set the
maximum number of contacts at a location: 
Here is an example where we set the duration of the time in the InfectedSymptoms state to the InfectedSevere state to 4 days:

.. code-block:: cpp

   model.parameters.get<mio::abm::TimeInfectedSymptomsToSevere>() = 4.;

We can also set the contact rates for specific age groups at a location:

.. code-block:: cpp

   model.get_location(work)
       .get_infection_parameters()
       .get<mio::abm::ContactRates>()[{age_group_15_to_34, age_group_15_to_34}] = 10.0;

For a full list of parameters, see `here <https://memilio.readthedocs.io/en/latest/api/file__home_docs_checkouts_readthedocs.org_user_builds_memilio_checkouts_latest_cpp_models_abm_parameters.h.html>`_.

Locations and persons
~~~~~~~~~~~~~~~~~~~~~

To add a location to the model, we have to specify the kind of location:

.. code-block:: cpp

   auto home = model.add_location(mio::abm::LocationType::Home);

People are added with an age. Then we have to assign them, so the model knows they can travel to this location:

.. code-block:: cpp

   auto person = model.add_person(home, age_group_0_to_4);
   person.set_assigned_location(home);

Note that adding the person to the model in one location does not mean that this location is in the list of assigned locations the person can visit afterwards.

For more complex location configurations, the model allows setting location-specific parameters:

.. code-block:: cpp

   // Add one social event with 5 maximum contacts (local)
   auto event = model.add_location(mio::abm::LocationType::SocialEvent);
   model.get_location(event).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
   
   // Increase aerosol transmission for all locations (global)
   model.parameters.get<mio::abm::AerosolTransmissionRates>() = 10.0;
   
   // Increase contact rate for specific age groups at a specific work location (local)
   auto work = model.add_location(mio::abm::LocationType::Work);
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

In this example, children are created in the age groups 0-4 and 5-14, while parents are created in the age groups 15-34 and 35-59, with equal weights respectively.

Testing strategies
~~~~~~~~~~~~~~~~~~

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
~~~~~~~~~~~~~~~~~~~~~~~

For infections to happen during the simulation, we have to initialize people with infections. Here, we iterate over all persons of the model and initialize them with random infection states according to a discrete distribution, i.e., 50% of persons are initialized as Susceptible, 30% as Exposed, etc.

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
~~~~~~~~~~~~~~~~~~~~~~

Here, we run the simulation:

.. code-block:: cpp

   auto t0   = mio::abm::TimePoint(0);
   auto tmax = t0 + mio::abm::days(30);
   auto sim  = mio::abm::Simulation(t0, std::move(model));
   
   // Simple simulation without data collection
   sim.advance(tmax);

Alternatively, if we want to track things in the simulation, we need to set up a
`history <https://github.com/SciCompMod/memilio/blob/main/cpp/memilio/io/README.md#the-history-object>`_, for example, to track all the Infection states of each simulation step into a Timeseries.

.. code-block:: cpp

   mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> history{
       Eigen::Index(mio::abm::InfectionState::Count)};

Then we can run the simulation with the history object and access the data through ``get_log()``:

.. code-block:: cpp

   sim.advance(tmax, history);
   auto log = history.get_log();

Finally, for example, we can print the data to a text file:

.. code-block:: cpp

   std::ofstream outfile("abm_minimal.txt");
   std::get<0>(log).print_table({"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4, outfile);
   std::cout << "Results written to abm_minimal.txt" << std::endl;

