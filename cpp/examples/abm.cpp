#include "epidemiology/abm/abm.h"

int main()
{
    auto world = epi::World();
    auto& home = world.add_node(epi::NodeType::Home);
    auto& school = world.add_node(epi::NodeType::School);
    auto& work = world.add_node(epi::NodeType::Work);
    auto& child1 = world.add_person(home, epi::InfectionState::Susceptible);
    auto& child2 = world.add_person(home, epi::InfectionState::Susceptible);
    auto& parent1 = world.add_person(home, epi::InfectionState::Exposed);
    auto& parent2 = world.add_person(home, epi::InfectionState::Susceptible);

    auto sim = epi::AbmSimulation(0, world);

    sim.advance(100);
}