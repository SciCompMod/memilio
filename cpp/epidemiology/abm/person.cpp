#include "epidemiology/abm/world.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/node.h"

namespace epi
{

Person::Person(Node& node, InfectionState state)
    : m_node(node)
    , m_state(state)
{
}

void Person::interact(double dt, const GlobalInfectionParameters& global_infection_params)
{
    auto state = m_state;
    auto new_state = state;

    if (state == InfectionState::Exposed) {
        if (m_time_until_carrier <= 0) {
            new_state = InfectionState::Carrier;
        }
        m_time_until_carrier -= dt;
    }
    else {
        new_state = m_node.get().interact(*this, dt);
        if (new_state == InfectionState::Exposed) {
            m_time_until_carrier = global_infection_params.incubation_time;
        }
    }

    m_state = new_state;
    if (state != new_state) {
        m_node.get().changed_state(*this, state);
    }
}

void Person::migrate_to(Node& node)
{
    m_node.get().remove_person(*this);
    m_node = node;
    m_node.get().add_person(*this);
}
} // namespace epi