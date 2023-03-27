#include "test_abm.h"

void add_infection_simple(mio::abm::Person& p, mio::abm::TimePoint t, mio::abm::InfectionState infection_state)
{
    mio::abm::GlobalInfectionParameters params;
    auto infection = mio::abm::Infection(static_cast<mio::abm::VirusVariant>(0), static_cast<mio::abm::AgeGroup>(0),
                                         params, t, infection_state);
    p.add_new_infection(mio::abm::Infection(static_cast<mio::abm::VirusVariant>(0), static_cast<mio::abm::AgeGroup>(0),
                                            params, t, infection_state));
}

mio::abm::Person create_person_simple(mio::abm::Location& location, mio::abm::InfectionState infection_state,
                                      mio::abm::AgeGroup age_group)
{
    mio::abm::Person p = mio::abm::Person(location, age_group);
    if (infection_state != mio::abm::InfectionState::Susceptible) {
        add_infection_simple(p, mio::abm::TimePoint(0), infection_state);
    }
    return p;
}

mio::abm::Person& add_person_simple(mio::abm::World& world, mio::abm::LocationId loc_id,
                                    mio::abm::InfectionState infection_state, mio::abm::AgeGroup age)
{
    mio::abm::Person& p = world.add_person(loc_id, age);
    if (infection_state != mio::abm::InfectionState::Susceptible) {
        add_infection_simple(p, mio::abm::TimePoint(0), infection_state);
    }
    return p;
}