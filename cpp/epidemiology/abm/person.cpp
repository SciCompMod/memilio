#include "epidemiology/abm/person.h"
#include "epidemiology/abm/node.h"

namespace epi
{

Person::Person(Node& node, InfectionState state)
    : m_node(node)
    , m_state(state)
{
}

void Person::interact(double dt)
{
    m_state = get_node().next_infection_state(*this, dt);
}

void Person::migrate_to(Node& node)
{
    get_node().remove_person(*this);
    m_node = node;
    get_node().add_person(*this);
}
} // namespace epi