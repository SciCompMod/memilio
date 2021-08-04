#include "epidemiology/abm/person.h"
#include "epidemiology/abm/world.h"
#include "epidemiology/abm/location.h"
#include "epidemiology/utils/random_number_generator.h"

namespace epi
{

Person::Person(LocationId id, InfectionState state, AbmAgeGroup age)
    : m_location_id(id)
    , m_assigned_locations((uint32_t)LocationType::Count, INVALID_LOCATION_INDEX)
    , m_state(state)
    , m_quarantine(false)
    , m_age(age)
    , m_time_at_location(std::numeric_limits<int>::max())
{
    m_random_workgroup = UniformDistribution<double>::get_instance()();
    m_random_schoolgroup = UniformDistribution<double>::get_instance()();
    if (state == InfectionState::Infected_Detected) {
        m_quarantine = true;
    }
}

Person::Person(Location& location, InfectionState state, AbmAgeGroup age)
    : Person({location.get_index(), location.get_type()}, state, age)
{
}

void Person::interact(TimeSpan dt, const GlobalInfectionParameters& global_infection_params, Location& loc)
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
        new_state = loc.interact(*this, dt, global_infection_params);
        if (new_state == InfectionState::Exposed) {
            m_time_until_carrier = hours(int(global_infection_params.get<IncubationPeriod>()[{this->m_age}] * 24));
        }
    }

    if (new_state == InfectionState::Infected_Detected) {
        m_quarantine = true;
    }
    else {
        m_quarantine = false;
    }

    m_state = new_state;
    if (state != new_state) {
        loc.changed_state(*this, state);
    }

    m_time_at_location += dt;
}

void Person::migrate_to(Location& loc_old, Location& loc_new)
{
    if (&loc_old!= &loc_new) {
        loc_old.remove_person(*this);
        m_location_id = {loc_new.get_index(), loc_new.get_type()};
        loc_new.add_person(*this);
        m_time_at_location = epi::TimeSpan(0);
    }
}

void Person::set_assigned_location(Location& location)
{
    m_assigned_locations[(uint32_t)location.get_type()] = location.get_index();
}

void Person::set_assigned_location (LocationId id)
{
    m_assigned_locations[(uint32_t)id.type] = id.index;
}

uint32_t Person::get_assigned_location_index (LocationType type) const
{
    return m_assigned_locations[(uint32_t)type];
}

bool Person::goes_to_work(TimePoint t, const AbmMigrationParameters& params) const
{
    return m_random_workgroup < params.get<WorkRatio>().get_matrix_at(t.days())[0];
}

bool Person::goes_to_school(TimePoint t, const AbmMigrationParameters& params) const
{
    return m_random_schoolgroup < params.get<SchoolRatio>().get_matrix_at(t.days())[0];
}
} // namespace epi
