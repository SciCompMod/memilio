/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#include "abm/time.h"
#include "memilio/utils/random_number_generator.h"
#include <algorithm>
#include <array>
#include <numeric>

namespace mio
{
namespace abm
{

/**
 * @brief Select a random transition from a list of possible transitions from the current state to others.
 * Each transition is represented by the new state and the probability of the transition, e.g. 
 * a pair {1, 0.5} is the transition to state 1 with rate 0.5.
 * Transition rates are not probabilities but the parameters of an exponential distribution.
 * One of the transitions happens if x < dt, where x is a sample from the exponential distribution Exp(S), 
 * S begin the sum of all rates. Which transition happens is determined by sampling from a discrete distribution
 * with the rates as weights. It's also possible that no transition happens in this time step.
 * In this case the current state is returned.
 * @tparam T Type that represents the states.
 * @tparam NumTransitions Number of possible transitions.
 * @param[in] current_state Current state before transition.
 * @param[in] dt Length of the time step.
 * @param[in] transitions Array of pairs of new states and their rates (probabilities).
 * @return New state from the list if transition happens, current_state otherwise.
 */
template <class T, size_t NumTransitions>
T random_transition(T current_state, TimeSpan dt, const std::pair<T, double> (&transitions)[NumTransitions])
{
    assert(std::all_of(std::begin(transitions), std::end(transitions),
                       [](auto& p) {
                           return p.second >= 0.0;
                       }) &&
           "transition rates must be non-negative");

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

} // namespace abm
} // namespace mio
