# Agent-Based Model

This module models and simulates the epidemic using an agent-based model (*ABM*) approach. Unlike the [SECIR](../ode_secir/README.md) compartmental model that uses a system of ODEs, this model simulates the spread of COVID-19 in a population with discrete persons (the agents) moving throughout locations in the world and interacting with (infecting) each other.

## Structure

The model consists of the following major classes:

1. Person: represents an agent of the model. A person has an ID, i.e. a unique number, an age, a location and a list with their assigned locations, i.e. the locations they visit during the simulation. They can perform tests and wear masks. Every person has lists with past and current infections and vaccinations.
2. Infection: collection of all information about a persons' infection, i.e. infectiousness, infection course, virus variant. The infection course is drawn stochastically from the infection states that are similar to the compartments of the SECIR model.
3. Location: represents places in the world where people meet and interact, e.g. home, school, work, social event sites. A location can be split into cells to model parts of a location, like classrooms in a school. Some infection parameters are location-specific and one can activate NPIs like mandatory masks or tests to enter the location.
4. World: collection of all persons and locations. It also holds information about the testing strategy of the simulation and holds the rules for the migration phase.
5. Simulation: runs the simulation and stores results.

## Simulation

The simulation runs in discrete time steps. Each step has two phases, an interaction phase and a migration phase.

### Interaction Phase

In this phase, each person interacts with the other persons at the same location. This interaction determines the transmission of the disease. A susceptible person can become infected by contact with an infected person. The probability of infection depends on a multitude of factors, such as the viral load and infectiousness of the infected and the immunity level of the susceptible person.

### Migration Phase

During the migration phase, each person may change their location. Migration follows complex [rules](../abm/migration_rules.cpp), considering the current location, time of day, and properties of the person (e.g. age). Some location changes are deterministic and regular (e.g. going to work), others are random (e.g. going to shopping or to a social event in the evening/on the weekend). When agents are infected, they are quarantined and cannot migrate. You can restrict some migration rules by allowing only a proportion of people to enter specific locations.

Another way of migration we use in the simulation of Braunschweig (simulations/abm_braunschweig.cpp) is using trips. A trip consists of the ID of the person that performs this trip, a time point when this trip is performed and where the person is heading to. At the beginning of the simulation, a list with all trips is initialized and followed during the simulation. There can be different trips on the weekend than during the week, but other than that, the agents do the same trips every day. As before, agents that are in quarantine or in the hospital cannot migrate.

## Get Started

This section gives an introduction to how to use the ABM and set up your own simulation. For a quick overview, can find a full example in the [ABM minimal example](../../examples/abm_minimal.cpp) and a more detailed Doxygen documentation [here](https://scicompmod.github.io/memilio/documentation/index.html ). For a guide on installation and running the simulations and examples, see this [README](../../README.md).

Every person in the ABM belongs to an AgeGroup, which we can define as the following. Note that every age group has to have values strictly smaller than num_age_groups.

```cpp
size_t num_age_groups         = 4;
const auto age_group_0_to_4   = mio::AgeGroup(0);
```

The initial empty world is created with the number of age groups:

```cpp
auto world = mio::abm::World(num_age_groups);
```

We can set plenty of general parameters, which you can find [here](../abm/parameters.h). Here is an example where we set the duration of the incubation period to 4 days:

```cpp
world.parameters.get<mio::abm::IncubationPeriod>() = 4.
```

To add a location to the world, we have to specify the kind of location.

```cpp
auto home = world.add_location(mio::abm::LocationType::Home);
```

People are added with an age. Then we have to assign them, so the world knows they can travel to this location.

```cpp
auto person = world.add_person(home, mio::AgeGroup(0));
person.set_assigned_location(home);
```

For adding more people to the world, we create households. A Household holds a vector with HouseholdMembers, i.e. a vector with weighted age distribution from which the age of the Persons belonging to this Household can be calculated. A Household and the number of times it exists is gathered in a Household Group.
For example, we have children who either belong to AgeGroup(0) or AgeGroup(1) with probability 0.5 in each case and parents which belong to AgeGroup(2) or AgeGroup(3) similarly. We then add households with a parent and a child and households with two parents and one child.

```cpp
auto child = mio::abm::HouseholdMember(num_age_groups);
child.set_age_weight(age_group_0_to_4, 1);
child.set_age_weight(age_group_0_to_4, 1);

auto parent = mio::abm::HouseholdMember(num_age_groups);
parent.set_age_weight(AGE_GROUP_15_TO_34, 1);
parent.set_age_weight(AGE_GROUP_35_TO_59, 1);

// Two-person household with one parent and one child.
auto twoPersonHousehold_group = mio::abm::HouseholdGroup();
auto twoPersonHousehold_full  = mio::abm::Household();
twoPersonHousehold_full.add_members(child, 1);
twoPersonHousehold_full.add_members(parent, 1);
twoPersonHousehold_group.add_households(twoPersonHousehold_full, n_households);
add_household_group_to_world(world, twoPersonHousehold_group);

```

During the simulation, people can get tested, and we have to specify the scheme for that:

```cpp
auto testing_min_time = mio::abm::days(1);
auto probability      = 0.5;
auto start_date       = mio::abm::TimePoint(0);
auto end_date         = mio::abm::TimePoint(0) + mio::abm::days(30);
auto test_type        = mio::abm::AntigenTest();
auto test_at_work     = std::vector<mio::abm::LocationType>{mio::abm::LocationType::Work};
auto testing_criteria_work =
    std::vector<mio::abm::TestingCriteria>{mio::abm::TestingCriteria({}, test_at_work, {})};
auto testing_scheme_work =
    mio::abm::TestingScheme(testing_criteria_work, testing_min_time, start_date, end_date, test_type, probability);
world.get_testing_strategy().add_testing_scheme(testing_scheme_work);
```

For some infections to happen during the simulation, we have to initialize people with infections.

```cpp
person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(), world.parameters, start_date, infection_state));
```

We can add restrictions for people after a specific date. For example, only 10% of the people go to social events after day 10.

```cpp
auto t_lockdown = mio::abm::TimePoint(0) + mio::abm::days(10);
mio::abm::close_social_events(t_lockdown, 0.9, world.parameters);
```

Then we run the simulation.

```cpp
sim.advance(mio::abm::TimePoint(0) + mio::abm::days(30));
```

## Current Limitations

Currently, a few things are not yet implemented, such as:

- Different trips for each day.
- Test and Trace functionality
