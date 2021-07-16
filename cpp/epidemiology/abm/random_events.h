#include "epidemiology/abm/time.h"
#include "epidemiology/utils/random_number_generator.h"
#include <algorithm>
#include <array>
#include <numeric>

namespace epi
{

/**
 * select a random transition from a list of possible transitions from the current state to others.
 * Each transition is represented by the new state and the probability of the transition, e.g. 
 * a pair {1, 0.5} is the transition to state 1 with rate 0.5.
 * Transition rates are not probabilities but the parameters of an exponential distribution.
 * One of the transitions happens if x < dt, where x is a sample from the exponential distribution Exp(S), 
 * S begin the sum of all rates. Which transition happens is determined by sampling from a discrete distribution
 * with the rates as weights. It's also possible that no transition happens in this time step.
 * In this case the current state is returned.
 * @tparam T type that represents the states
 * @tparam NumTransitions number of possible transitions
 * @param current_state current state before transitions
 * @param dt length of the time step
 * @param transitions array of pairs of new states and their rates (probabilities)
 * @return new state from the list if transition happens, current_state otherwise
 */
template <class T, size_t NumTransitions>
T random_transition(T current_state, TimeSpan dt, const std::pair<T, double> (&transitions)[NumTransitions])
{
    assert(std::all_of(std::begin(transitions), std::end(transitions), [](auto& p) {
        return p.second >= 0.0;
    }) && "transition rates must be non-negative");

    //check if any transition happens using exponential distribution with the sum of all transition rates
    auto sum = std::accumulate(std::begin(transitions), std::end(transitions), 0.0, [](auto&& a, auto&& t) {
        return a + t.second;
    });
    if (sum <= 0) { //no transitions or all transitions have rate zero
        return current_state;
    }
    auto v = ExponentialDistribution<double>::get_instance()(sum);
    if (v < dt.days()) {
        //pick one of the possible transitions using discrete distribution
        std::array<double, NumTransitions> rates;
        std::transform(std::begin(transitions), std::end(transitions), rates.begin(), [](auto&& t) {
            return t.second;
        });
        auto random_idx = DiscreteDistribution<size_t>::get_instance()(rates);
        return transitions[random_idx].first;
    }

    return current_state;
}

} // namespace epi