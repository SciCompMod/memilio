#include "epidemiology/abm/node.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/rng.h"

namespace epi
{

Node::Node(NodeType type)
    : m_type(type)
    , m_subpopulations{}
{
}

namespace
{
    template <int NumTransitions>
    InfectionState random_transition(InfectionState old_state, double dt,
                                     const std::pair<InfectionState, double> (&transition_rates)[NumTransitions])
    {
        std::pair<InfectionState, double> transition_rates_partial_sum[NumTransitions];
        std::partial_sum(std::begin(transition_rates), std::end(transition_rates),
                         std::begin(transition_rates_partial_sum), [](auto&& rate_sum, auto&& rate) {
                             auto tmp = rate;
                             tmp.second += rate_sum.second;
                             return tmp;
                         });
        auto sum = transition_rates_partial_sum[NumTransitions - 1].second;

        if (sum <= 0)
        {
            return old_state;
        }

        auto uni_dist = std::uniform_real_distribution<double>(0, 1);
        auto v = -std::log(uni_dist(thread_local_rng())) / sum; //determines if there is a transition in this step
        auto u = uni_dist(thread_local_rng()); //determines which transition happens

        if (v < dt) {
            auto transition =
                *std::lower_bound(std::begin(transition_rates_partial_sum), std::end(transition_rates_partial_sum),
                                  u * sum, [](auto&& rate, auto&& u) {
                                      return rate.second < u;
                                  }); //scaling u makes sure that lower bound on the partial sums succeeds
            return transition.first;
        }

        return old_state;
    }
} // namespace

InfectionState Node::interact(const Person& person, double dt) const
{
    auto state = person.get_infection_state();
    switch (state) {
    case InfectionState::Susceptible:
        return random_transition(state, dt, {{InfectionState::Exposed, dt * m_cached_exposure_rate}});        
    case InfectionState::Carrier:
        return random_transition(state, dt,
            {{InfectionState::Infected_Detected, m_parameters.detect_infection * m_parameters.carrier_to_infected},
             {InfectionState::Infected_Undetected, (1 - m_parameters.detect_infection) * m_parameters.carrier_to_infected},
             {InfectionState::Recovered_Carrier, m_parameters.carrier_to_recovered}});             
    case InfectionState::Infected_Detected: //fallthrough!
    case InfectionState::Infected_Undetected:
        return random_transition(state, dt,
            {{InfectionState::Recovered_Infected, m_parameters.infected_to_recovered},
             {InfectionState::Dead, m_parameters.infected_to_dead * m_parameters.death_factor}});
    case InfectionState::Recovered_Carrier: //fallthrough!
    case InfectionState::Recovered_Infected:
        return random_transition(state, dt, {{InfectionState::Susceptible, m_parameters.recovered_to_susceptible}});
    default:
        return state; //some states don't transition
    }
}

void Node::end_step(double dt)
{
    //cache for next step so it stays constant during the step while subpopulations change
    //otherwise we would have to cache all state changes during a step which uses more memory
    auto num_carriers = get_subpopulation(InfectionState::Carrier);
    auto num_infected =
        get_subpopulation(InfectionState::Infected_Detected) + get_subpopulation(InfectionState::Infected_Undetected);
    m_cached_exposure_rate = std::min(m_parameters.effective_contacts, double(m_num_persons)) / m_num_persons *
                             (m_parameters.susceptible_to_carrier_by_carrier * num_carriers +
                              m_parameters.susceptible_to_carrier_by_infected * num_infected);
}

void Node::add_person(Person& p)
{
    ++m_num_persons;
    InfectionState s = p.get_infection_state();
    change_subpopulation(s, +1);
}

void Node::remove_person(Person& p)
{
    --m_num_persons;
    InfectionState s = p.get_infection_state();
    change_subpopulation(s, -1);
}

void Node::changed_state(const Person& p, InfectionState old_state)
{
    change_subpopulation(old_state, -1);
    change_subpopulation(p.get_infection_state(), +1);
}

void Node::change_subpopulation(InfectionState s, int delta)
{
    m_subpopulations[size_t(s)] += delta;
}

int Node::get_subpopulation(InfectionState s) const
{
    return m_subpopulations[size_t(s)];
}
} // namespace epi
