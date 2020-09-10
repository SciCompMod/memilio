# Agent Based Model

This module models and simulates the epidemic using an agent based model (*ABM*) approach. Unlike the SECIR compartment model that uses an ODE, this model simulates the spread of COVID19 in a population with actual persons (the agents) moving throughout the world and interacting with each other.

## Structure

The model consists of the following major classes:
1. Person: represents the agents of the model. A person has an assigned infection state (similar to the compartments of the SECIR model) and location. 
2. Location: represents locations in the world where people meet and interact, e.g. home, school, work, social events.
3. World: collection of all persons and locations.
4. Simulation: run the simulation and store results.

## Simulation

The simulation runs in discrete time steps. Each step is in two phases, an interaction phase and a migration phase. First each person interacts with the other persons at the same location. This interaction determines their infection state transitions. A susceptible person may become exposed by contact with a carrier, an asymptomatic carrier can become symptomatic (infected), etc. The transitions are determined by chance. The probability of a transition depends on the parameters of the infection and the population at this location. During the migration phase, each person may change location. Migration follows complex rules, taking into account the current location, time of day, and properties of the person (e.g. age) (NOT YET IMPLEMENTED). The result of the simulation is for each time step the count of persons in each infection state at that time.

## Example
An example can also be found at examples/abm.cpp
```
    //setup the infection parameters
    epi::GlobalInfectionParameters params;
    params.incubation_period = 4.7;
    //...

    //create the world with some nodes and persons
    epi::World world = epi::World(params);
    epi::Location& home = world.add_location(epi::LocationType::Home);
    epi::Location& school = world.add_location(epi::LocationType::School);
    epi::Location& work = world.add_location(epi::LocationType::Work);
    epi::Person& p1 = world.add_person(home, epi::InfectionState::Susceptible);
    epi::Person& p2 = world.add_person(home, epi::InfectionState::Susceptible);
    epi::Person& p3 = world.add_person(home, epi::InfectionState::Carrier);
    epi::Person& p4 = world.add_person(home, epi::InfectionState::Susceptible);

    //setup the simulation
    double t0   = 0;
    double tmax = 100;
    epi::AbmSimulation sim  = epi::AbmSimulation(t0, std::move(world));

    sim.advance(tmax);

    //handle result
    epi::TimeSeries<double>& result = sim.get_result();
    //...
```
