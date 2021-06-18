#include "epidemiology/abm/abm.h"

int main()
{
    auto world    = epi::World();
    auto home    = world.add_location(epi::LocationType::Home);
    auto school  = world.add_location(epi::LocationType::School);
    auto work    = world.add_location(epi::LocationType::Work);
    auto child1  = world.add_person(home, epi::InfectionState::Susceptible);
    auto child2  = world.add_person(home, epi::InfectionState::Susceptible);
    auto parent1 = world.add_person(home, epi::InfectionState::Carrier);
    auto parent2 = world.add_person(home, epi::InfectionState::Susceptible);

    auto t0   = 0;
    auto tmax = 100;
    auto sim  = epi::AbmSimulation(epi::TimePoint(t0), std::move(world));

    sim.advance(epi::TimePoint(tmax));

    std::cout << "Ran ABM from " << t0 << " to " << 100 << '\n';
    for (auto&& v : sim.get_result()) {
        std::cout << v.transpose() << '\n';
    }
    std::cout << std::endl;
}
