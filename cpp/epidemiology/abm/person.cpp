#include "epidemiology/abm/person.h"
#include "epidemiology/abm/world.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/abm/random_number_generator.h"

namespace epi
{

Person::Person(Location& location, InfectionState state, AbmAgeGroup age)
    : m_location(location)
    , m_state(state)
    , m_age(age)
    , m_time_at_location(std::numeric_limits<int>::max())
{
}

void Person::interact(TimeSpan dt, const GlobalInfectionParameters& global_infection_params)
{
    auto state     = m_state;
    auto new_state = state;

    if (state == InfectionState::Exposed) {
        if (m_time_until_carrier <= TimeSpan(0)) {
            new_state = InfectionState::Carrier;
        }
        m_time_until_carrier -= dt;
    }
    else {
        new_state = m_location.get().interact(*this, dt, global_infection_params);
        if (new_state == InfectionState::Exposed) {
            m_time_until_carrier = hours(int(global_infection_params.get<IncubationPeriod>()[{this->m_age}] * 24));
        }
    }

    m_state = new_state;
    if (state != new_state) {
        m_location.get().changed_state(*this, state);
    }

    m_time_at_location += dt;
}

void Person::migrate_to(Location& location)
{
    if (&location != &m_location.get()) {
        m_location.get().remove_person(*this);
        m_location = location;
        m_location.get().add_person(*this);
        m_time_at_location = TimeSpan(0);
    }
}
} // namespace epi
