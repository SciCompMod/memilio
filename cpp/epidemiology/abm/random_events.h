#include "epidemiology/abm/random_number_generator.h"
#include "epidemiology/abm/time.h"
#include <algorithm>
#include <array>
#include <numeric>

namespace epi
{

/**
 * select a random transition from a list of possible transitions.
 * each transition is represented by the new state and the probability of the transition.
 * it's also possible that no transition happens in this time step.
 * in this case the current state is returned.
 * input rates don't need to sum up to probability 1, the function performs normalisation.
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
    auto sum = std::accumulate(std::begin(transitions), std::end(transitions), 0.0, [](auto&& a, auto&& t) {
        return a + t.second;
    });

    if (sum <= 0) {
        return current_state;
    }

    //time between transitions is exponentially distributed
    auto v = ExponentialDistribution<double>::get_instance()(sum);
    if (v < dt.days()) {
        //pick a random transition
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