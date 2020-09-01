#include "epidemiology/abm/location.h"
#include "epidemiology/abm/person.h"
#include "epidemiology/abm/random_number_generator.h"

namespace epi
{

Location::Location(LocationType type)
    : m_type(type)
    , m_subpopulations{}
{
}

namespace
{
    /**
     * select a random transition from a list of possible transitions.
     * each transition is represented by the new state and the probability of the transition.
     * it's also possible that no transition happens in this time step.
     * in this case the old state is returned.
     * input rates don't need to sum up to probability 1, the function performs normalisation.
     * @param old_state previous state
     * @param dt length of the time step
     * @param transition_rates array of pairs of new states and their rates (probabilities)
     * @return new state if transition happens, old_state otherwise
     */
    template <int NumTransitions>
    InfectionState random_transition(InfectionState old_state, double dt,
                                     const std::pair<InfectionState, double> (&transition_rates)[NumTransitions])
    {
        std::pair<InfectionState, double> transition_rates_partial_sum[NumTransitions];
        std::partial_sum(std::begin(transition_rates), std::end(transition_rates),
                         std::begin(transition_rates_partial_sum), [](auto&& rate_sum, auto&& rate) {
                             //sum up the rates, but keep the states as they are
                             auto tmp = rate;
                             tmp.second += rate_sum.second;
                             return tmp;
                         });
        auto sum = transition_rates_partial_sum[NumTransitions - 1].second;

        if (sum <= 0) {
            return old_state;
        }

        auto uni_dist = std::uniform_real_distribution<double>(0, 1);
        auto v        = -std::log(uni_dist(thread_local_rng())) /
                 sum; //determines if there is a transition in this step (poisson distributed?)
        auto u = uni_dist(thread_local_rng()); //determines which transition happens

        if (v < dt) {
            auto transition =
                *std::lower_bound(std::begin(transition_rates_partial_sum), std::end(transition_rates_partial_sum),
                                  u * sum, [](auto&& rate, auto&& u) {
                                      return rate.second < u;
                                  }); //scaling u to map the transition rates onto uniform distribution in [0, 1)
            return transition.first;
        }

        return old_state;
    }
} // namespace

InfectionState Location::interact(const Person& person, double dt, const GlobalInfectionParameters& global_params) const
{
    auto state = person.get_infection_state();
    switch (state) {
    case InfectionState::Susceptible:
        return random_transition(state, dt, {{InfectionState::Exposed, dt * m_cached_exposure_rate}});
    case InfectionState::Carrier:
        return random_transition(
            state, dt,
            {{InfectionState::Infected_Detected, global_params.detect_infection * global_params.carrier_to_infected},
             {InfectionState::Infected_Undetected,
              (1 - global_params.detect_infection) * global_params.carrier_to_infected},
             {InfectionState::Recovered_Carrier, global_params.carrier_to_recovered}});
    case InfectionState::Infected_Detected: //fallthrough!
    case InfectionState::Infected_Undetected:
        return random_transition(state, dt,
                                 {{InfectionState::Recovered_Infected, global_params.infected_to_recovered},
                                  {InfectionState::Dead, global_params.infected_to_dead * m_parameters.death_factor}});
    case InfectionState::Recovered_Carrier: //fallthrough!
    case InfectionState::Recovered_Infected:
        return random_transition(state, dt, {{InfectionState::Susceptible, global_params.recovered_to_susceptible}});
    default:
        return state; //some states don't transition
    }
}

void Location::begin_step(double dt, const GlobalInfectionParameters& global_params)
{
    //cache for next step so it stays constant during the step while subpopulations change
    //otherwise we would have to cache all state changes during a step which uses more memory
    auto num_carriers = get_subpopulation(InfectionState::Carrier);
    auto num_infected =
        get_subpopulation(InfectionState::Infected_Detected) + get_subpopulation(InfectionState::Infected_Undetected);
    m_cached_exposure_rate = std::min(m_parameters.effective_contacts, double(m_num_persons)) / m_num_persons *
                             (global_params.susceptible_to_exposed_by_carrier * num_carriers +
                              global_params.susceptible_to_exposed_by_infected * num_infected);
}

void Location::add_person(const Person& p)
{
    ++m_num_persons;
    InfectionState s = p.get_infection_state();
    change_subpopulation(s, +1);
}

void Location::remove_person(const Person& p)
{
    --m_num_persons;
    InfectionState s = p.get_infection_state();
    change_subpopulation(s, -1);
}

void Location::changed_state(const Person& p, InfectionState old_state)
{
    change_subpopulation(old_state, -1);
    change_subpopulation(p.get_infection_state(), +1);
}

void Location::change_subpopulation(InfectionState s, int delta)
{
    m_subpopulations[size_t(s)] += delta;
}

int Location::get_subpopulation(InfectionState s) const
{
    return m_subpopulations[size_t(s)];
}

Eigen::Ref<const Eigen::VectorXi> Location::get_subpopulations() const
{
    return Eigen::Map<const Eigen::VectorXi>(m_subpopulations.data(), m_subpopulations.size());
}
} // namespace epi
