# Agent Based Model

This module models and simulates the epidemic using an agent based model (*ABM*) approach. Unlike the SECIR compartment model that uses an ODE, this model simulates the spread of COVID19 in a population with actual persons (the agents) moving throughout the world and interacting with each other.

## Get Started

The model consists of a world with persons and locations. Every person belongs to an AgeGroup which we can define.
```
    size_t NUM_AGE_GROUPS         = 4;
    const auto AGE_GROUP_0_TO_4   = mio::AgeGroup(NUM_AGE_GROUPS - 4);
    const auto AGE_GROUP_5_TO_14  = mio::AgeGroup(NUM_AGE_GROUPS - 3);
    const auto AGE_GROUP_15_TO_34 = mio::AgeGroup(NUM_AGE_GROUPS - 2);
    const auto AGE_GROUP_35_TO_59 = mio::AgeGroup(NUM_AGE_GROUPS - 1);

    // Create the world with 4 age groups.
    auto world = mio::abm::World(NUM_AGE_GROUPS);
```
We set some parameters and check if they satisfy their constraints:

```
    world.parameters.get<mio::abm::IncubationPeriod>() = 4.;

    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    world.parameters.get<mio::abm::AgeGroupGotoSchool>() = {AGE_GROUP_5_TO_14};
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 and 35-59)
    world.parameters.get<mio::abm::AgeGroupGotoWork>() = {AGE_GROUP_15_TO_34, AGE_GROUP_35_TO_59};

    // Check if the parameters satisfy their contraints.
    world.parameters.check_constraints();
```
For adding people to the world we create households. A Household holds a vector with HouseholdMembers, i.e. a vector with weighted age distribution from which the age of the Persons belonging to this Household can be calculated. A Household and the amount of times it exists is gathered in a HouseholdGroup.
For example, we have children which either belong to AgeGroup(0) or AgeGroup(1) with probability 0.5 in each case and parents which belong to AgeGroup(2) or AgeGroup(3) similarly. We then add households with a parent and a child and households with two parents and one child. 
```
    int n_households = 3;

    auto child = mio::abm::HouseholdMember(NUM_AGE_GROUPS);
    child.set_age_weight(AGE_GROUP_0_TO_4, 1);
    child.set_age_weight(AGE_GROUP_0_TO_4, 1);

    auto parent = mio::abm::HouseholdMember(NUM_AGE_GROUPS);
    parent.set_age_weight(AGE_GROUP_15_TO_34, 1);
    parent.set_age_weight(AGE_GROUP_35_TO_59, 1);

    // Two-person household with one parent and one child.
    auto twoPersonHousehold_group = mio::abm::HouseholdGroup();
    auto twoPersonHousehold_full  = mio::abm::Household();
    twoPersonHousehold_full.add_members(child, 1);
    twoPersonHousehold_full.add_members(parent, 1);
    twoPersonHousehold_group.add_households(twoPersonHousehold_full, n_households);
    add_household_group_to_world(world, twoPersonHousehold_group);

    // Three-person household with two parent and one child.
    auto threePersonHousehold_group = mio::abm::HouseholdGroup();
    auto threePersonHousehold_full  = mio::abm::Household();
    threePersonHousehold_full.add_members(child, 1);
    threePersonHousehold_full.add_members(parent, 2);
    threePersonHousehold_group.add_households(threePersonHousehold_full, n_households);
    add_household_group_to_world(world, threePersonHousehold_group);
```
Now we already have the persons of our world so we add the locations. For each location we set the maximum number of people that a person can infect while being at this location.
```
    // Add one social event with 5 maximum contacts.
    auto event = world.add_location(mio::abm::LocationType::SocialEvent);
    world.get_individualized_location(event).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    // Add hospital and ICU with 5 maximum contacs.
    auto hospital = world.add_location(mio::abm::LocationType::Hospital);
    world.get_individualized_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto icu = world.add_location(mio::abm::LocationType::ICU);
    world.get_individualized_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    // Add one supermarket, maximum constacts are assumed to be 20.
    auto shop = world.add_location(mio::abm::LocationType::BasicsShop);
    world.get_individualized_location(shop).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    // At every school, the maximum contacts are 20.
    auto school = world.add_location(mio::abm::LocationType::School);
    world.get_individualized_location(school).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    // At every workplace, maximum contacts are 10.
    auto work = world.add_location(mio::abm::LocationType::Work);
    world.get_individualized_location(work).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
```
During the simulation people can get tested at work (and do this with probability 0.5) until day 30. 
```
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
For some infections to happen during the simulation we have to initialize people with infections. Therefore we chose random infection states and add an infection if a person is not susceptible.
```
    auto persons = world.get_persons();
    for (auto& person : persons) {
        mio::abm::InfectionState infection_state =
            (mio::abm::InfectionState)(rand() % ((uint32_t)mio::abm::InfectionState::Count - 1));
        auto rng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         world.parameters, start_date, infection_state));
        }
    }
```
Then we assign the persons to the locations above so that the persons always use the same location, i.e. a child always goes to the same school.
<!--- Soll das nicht lieber oben direkt nach der Erstellung der Locations? --->
```
    for (auto& person : persons) {
        //assign shop and event
        person.set_assigned_location(event);
        person.set_assigned_location(shop);
        //assign hospital and ICU
        person.set_assigned_location(hospital);
        person.set_assigned_location(icu);
        //assign work/school to people depending on their age
        if (person.get_age() == AGE_GROUP_0_TO_4) {
            person.set_assigned_location(school);
        }
        if (person.get_age() == AGE_GROUP_15_TO_34 || person.get_age() == AGE_GROUP_35_TO_59) {
            person.set_assigned_location(work);
        }
    }
```
We can add restrictions for people after a specific date. In this example only 10% of the people go to social events after day 10.
```
    auto t_lockdown = mio::abm::TimePoint(0) + mio::abm::days(10);
    mio::abm::close_social_events(t_lockdown, 0.9, world.parameters);
```
Then we set the duration of the simulation. We run the simulation and save the results.
```
    auto t0   = mio::abm::TimePoint(0);
    auto tmax = mio::abm::TimePoint(0) + mio::abm::days(30);
    auto sim  = mio::abm::Simulation(t0, std::move(world));

    sim.advance(tmax);

    write_results_to_file(sim)
```

## Structure

The model consists of the following major classes:
1. Person: represents the agents of the model. A person has an ID, i.e. a unique number, an age, a location and a list with their assigned locations, i.e. the locations they visit during the simulation. Every person has lists with experienced infections and vaccinations. They can perform tests and wear masks.
2. Infection: collection of all information about a persons infection, i.e. infectiousness, infection course, virus variant. The infection course is drawn stochastically from the infection states that are similar to the compartments of the SECIR model. 
3. Location: represents locations in the world where people meet and interact, e.g. home, school, work, social events. A location can be split into cells to model parts of a location like classrooms in a school. Some infection parameters are location specific and one can activate NPIs like mandatory masks or tests to enter the location.
4. World: collection of all persons and locations. It also holds information about the testing strategy of the simulation and holds the rules for the migration phase.
5. Simulation: run the simulation and store results.

## Simulation

The simulation runs in discrete time steps. Each step is in two phases, an interaction phase and a migration phase. 

First each person interacts with the other persons at the same location. This interaction determines the transmission of the desease. A susceptible person can become infected by contact with an infected person. The probability of infection depends on a multitude of factors, such as the viral load and infectiousness of the infectee and the immunity level of the susceptible person.

During the migration phase, each person may change location. Migration follows complex rules, taking into account the current location, time of day, and properties of the person (e.g. age). Some location changes are deterministic and regular (e.g. going to work), others are random (e.g. going to shopping or to a social event in the evening/on the weekend).

Another way of migration we use in the simulation of Braunschweig (simulations/abm_braunschweig.cpp) is using trips. A trip consists of the ID of the person that performs this trip, a timepoint when this trip is performed and where the person is heading to. In the beginning of the simulation a list with all trips is initialized and followed during the simulation. There can be different trips on the weekend than during the week but other than that the agents do the same trips every day assuming they are not in quarantine or in the hospital.

The result of the simulation is for each time step the count of persons in each infection state at that time.

## Example

An example can also be found at simulations/abm.cpp
```
    size_t NUM_AGE_GROUPS         = 4;
    const auto AGE_GROUP_0_TO_4   = mio::AgeGroup(NUM_AGE_GROUPS - 4);
    const auto AGE_GROUP_5_TO_14  = mio::AgeGroup(NUM_AGE_GROUPS - 3);
    const auto AGE_GROUP_15_TO_34 = mio::AgeGroup(NUM_AGE_GROUPS - 2);
    const auto AGE_GROUP_35_TO_59 = mio::AgeGroup(NUM_AGE_GROUPS - 1);

    // Create the world with 4 age groups.
    auto world = mio::abm::World(NUM_AGE_GROUPS);
    //setup the parameters
    world.parameters.get<mio::abm::IncubationPeriod>() = 4.;
    //...

    //create the world with some nodes and persons
    mio::abm::Location& home = world.add_location(mio::abm::LocationType::Home);
    mio::abm::Location& school = world.add_location(mio::abm::LocationType::School);
    mio::abm::Location& work = world.add_location(mio::abm::LocationType::Work);
    mio::abm::Person& p1 = world.add_person(home, mio::abm::InfectionState::Susceptible);
    mio::abm::Person& p2 = world.add_person(home, mio::abm::InfectionState::Susceptible);
    mio::abm::Person& p3 = world.add_person(home, mio::abm::InfectionState::InfectedNoSymptoms);
    mio::abm::Person& p4 = world.add_person(home, mio::abm::InfectionState::Susceptible);

    //setup the simulation
    double t0   = 0;
    double tmax = 100;
    mio::abm::Simulation sim  = mio::abm::Simulation(t0, std::move(world));

    sim.advance(tmax);

    //handle result
    mio::TimeSeries<double>& result = sim.get_result();
    //...
```
