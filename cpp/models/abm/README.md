# Agent-Based Model

This module models and simulates the epidemic using an agent-based model (*ABM*) approach. Unlike compartmental models like the [SECIR](../ode_secir/README.md) model that uses a system of ODEs, this model simulates the spread of COVID-19 in a population with discrete persons (the agents) moving throughout locations in the model and interacting with (infecting) each other.

## Structure

The model consists of the following major classes:

1. Person: represents an agent of the model. A person has an ID, i.e. a unique number, an age, a location and a list with their assigned locations, i.e. the locations it can visit during the simulation. They can perform tests and wear masks. Every person has lists with past and current infections and vaccinations.
2. Infection: collection of all information about a persons' infection, i.e. infectiousness, infection course, virus variant. The infection course is drawn stochastically from the infection states that are similar to the compartments of the SECIR model.
3. Location: represents places in the model where people meet and interact, e.g. home, school, work, social event sites. A location can be split into cells to model parts of a location, like classrooms in a school. Some infection parameters are location-specific and one can activate NPIs like mandatory masks or tests to enter the location.
4. Model: collection of all persons and locations. It also holds information about the testing strategy of the simulation and holds the rules for the mobility phase.
5. Simulation: runs the simulation and stores results.

## Simulation

The simulation runs in discrete time steps. Each step has two phases, an interaction phase and a mobility phase.

### Interaction Phase

In this phase, each person interacts with the other persons at the same location. This interaction determines the transmission of the disease. A susceptible person can become infected by contact with an infected person. The probability of infection depends on a multitude of factors, such as the viral load and infectiousness of the infected and the immunity level of the susceptible person at the time of transmission.

### Mobility Phase

During the mobility phase, each person may change their location. Mobility follows complex [rules](../abm/mobility_rules.cpp), considering the current location, time of day, and properties of the person (e.g. age). Some location changes are deterministic and regular (e.g. going to work), others are random (e.g. going to shopping or to a social event in the evening/on the weekend). When agents are tested positive, they are quarantined and cannot leave their home, unless their infection becomes worse and they have to go to the hospital or the ICU. In general, if an agent suffers from severe or critical symptoms, it will move to the hospital or ICU with highest priority. Some mobility rules can be restricted by allowing only a proportion of people to enter specific locations.

Another way of mobility we use in the simulation of Braunschweig (simulations/abm_braunschweig.cpp) is using trips. A trip consists of the ID of the person that performs this trip, a time point when this trip is performed and where the person is heading to. At the beginning of the simulation, a list with all trips is initialized and followed during the simulation. There can be different trips on the weekend than during the week, but other than that, the agents do the same trips every day. As before, agents that are in quarantine or in the hospital cannot change their location.

## Get Started

This section gives an introduction to how to use the ABM and set up your own simulation. For a quick overview, can find a full example in the [ABM minimal example](../../examples/abm_minimal.cpp) and a more detailed Doxygen documentation [here](https://scicompmod.github.io/memilio/documentation/index.html ). For a guide on installation and running the simulations and examples, see this [README](../../README.md).

Every person in the ABM belongs to an AgeGroup, which we can define as follows:  

```cpp  
size_t num_age_groups         = 4;  
const auto age_group_0_to_4   = mio::AgeGroup(0);  
const auto age_group_5_to_14  = mio::AgeGroup(1);  
...                           = ...  
```  

Note that every age group has to have values strictly smaller than the number of age groups `num_age_groups`.  
With this number we create an empty model:  

```cpp
auto model = mio::abm::Model(num_age_groups);
```

We can set several general parameters, which you can find [here](../abm/parameters.h). Here is an example where we set the duration of the incubation period to 4 days:

```cpp
model.parameters.get<mio::abm::IncubationPeriod>() = 4.
```

To add a location to the model, we have to specify the kind of location.

```cpp
auto home = model.add_location(mio::abm::LocationType::Home);
```

People are added with an age. Then we have to assign them, so the model knows they can travel to this location.

```cpp
auto person = model.add_person(home, age_group_0_to_4);
person.set_assigned_location(home);
```

For adding more people to the model, we create households. A Household holds a vector of HouseholdMembers, which in turn hold a weighted distribution, such that we can randomly draw the age of each Person belonging to the Household. To manage multiple Households of the same type, we can use a HouseholdGroup.
In our example, we categorize individuals into two groups: children and parents.

- Children: They can either belong to AgeGroup(0) or AgeGroup(1). The probability of a child belonging to either group is 0.5.
- Parents: They can either belong to AgeGroup(2) or AgeGroup(3). The probability of a parent belonging to either group is also 0.5.

We then form households in two ways:

1. Households with one parent and one child.
2. Households with two parents and one child.

```cpp
auto child = mio::abm::HouseholdMember(num_age_groups);
child.set_age_weight(age_group_0_to_4, 1);
child.set_age_weight(age_group_0_to_4, 1);

auto parent = mio::abm::HouseholdMember(num_age_groups);
parent.set_age_weight(age_groups_15_to_34, 1);
parent.set_age_weight(age_groups_35_to_59, 1);

// Two-person household with one parent and one child.
auto twoPersonHousehold_group = mio::abm::HouseholdGroup();
auto twoPersonHousehold_full  = mio::abm::Household();
twoPersonHousehold_full.add_members(child, 1);
twoPersonHousehold_full.add_members(parent, 1);
twoPersonHousehold_group.add_households(twoPersonHousehold_full, n_households);
add_household_group_to_model(model, twoPersonHousehold_group);

```

During the simulation, people can get tested, and we have to specify the scheme for that:

```cpp
auto probability      = 0.5;
auto start_date       = mio::abm::TimePoint(0);
auto end_date         = mio::abm::TimePoint(0) + mio::abm::days(30);
auto test_type        = mio::abm::AntigenTest();
auto test_at_work     = std::vector<mio::abm::LocationType>{mio::abm::LocationType::Work};
auto testing_criteria_work =
    std::vector<mio::abm::TestingCriteria>{mio::abm::TestingCriteria({}, test_at_work, {})};
auto testing_scheme_work =
    mio::abm::TestingScheme(testing_criteria_work, start_date, end_date, test_type, probability);
model.get_testing_strategy().add_testing_scheme(testing_scheme_work);
```

For some infections to happen during the simulation, we have to initialize people with infections.

```cpp
person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(), model.parameters, start_date, infection_state));
```

We can add restrictions for people after a specific date. For example, only 10% of the people go to social events after day 10.

```cpp
auto t_lockdown = mio::abm::TimePoint(0) + mio::abm::days(10);
mio::abm::close_social_events(t_lockdown, 0.9, model.parameters);
```

Then we run the simulation.

```cpp
sim.advance(mio::abm::TimePoint(0) + mio::abm::days(30));
```

Alternitavely, if we want to track things in the simulation, we need to set up a [history](../../memilio/io/README.md#the-history-object), for example, to track all the Infection states of each simulation step.

```cpp
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> history;
```

Then we can run the simulation with the history object and one can access the data it through get_log.

```cpp
    sim.advance(tmax, history);
    auto log = history.get_log();
```

Finally we can print the data to a text file.

```cpp
    std::ofstream outfile("abm_minimal.txt");
    std::get<0>(log).print_table({"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4, outfile);
    std::cout << "Results written to abm_minimal.txt" << std::endl;
```

## Current Limitations

Currently, a few things are not yet implemented, such as:

- Different trips for each day.
- Test and Trace functionality
