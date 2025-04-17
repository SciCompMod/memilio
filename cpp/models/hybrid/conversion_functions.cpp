/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Julia Bicker
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

#include "hybrid/conversion_functions.h"
#include "d_abm/single_well.h"
#include "hybrid/infection_state.h"
#include "memilio/geography/regions.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"
#include <algorithm>
#include <cstddef>
#include <numeric>

namespace mio
{
namespace hybrid
{
template <>
void convert_model(const dabm::Simulation<SingleWell<InfectionState>>& current_model,
                   smm::Simulation<1, InfectionState>& target_model)
{
    auto& current_result = current_model.get_result();
    auto& target_result  = target_model.get_result();
    assert(current_result.get_last_time() >= target_result.get_last_time() &&
           "Conversion from dabm to smm not possible because last smm time point is bigger than last dabm time point.");
    if (target_result.get_last_time() < current_result.get_last_time()) {
        target_result.add_time_point(current_result.get_last_time());
    }
    //update result timeseries
    target_result.get_last_value() = current_result.get_last_value();
    //update model populations
    for (int i = 0; i < (int)InfectionState::Count; ++i) {
        target_model.get_model().populations[{regions::Region(0), InfectionState(i)}] =
            current_result.get_last_value()[target_model.get_model().populations.get_flat_index(
                {regions::Region(0), InfectionState(i)})];
    }
}

template <>
void convert_model(const smm::Simulation<1, InfectionState>& current_model,
                   dabm::Simulation<SingleWell<InfectionState>>& target_model)
{
    auto& current_result = current_model.get_result();
    auto& target_result  = target_model.get_result();
    assert(current_result.get_last_time() >= target_result.get_last_time() &&
           "Conversion from dabm to smm not possible because last smm time point is bigger than last dabm time point.");
    if (target_result.get_last_time() < current_result.get_last_time()) {
        target_result.add_time_point(current_result.get_last_time());
    }
    //Update result timeseries
    target_result.get_last_value() = current_result.get_last_value();

    //Update agents' infection state and sample agents position
    auto current_pop = current_result.get_last_value().eval();
    double total_pop = std::accumulate(current_pop.begin(), current_pop.end(), 0.0);
    SWPositionSampler pos_rng{{-1, -1}, {1, 1}, 0.1};
    auto& state_rng = DiscreteDistribution<size_t>::get_instance();
    auto& abm_pop   = target_model.get_model().populations;
    if (abm_pop.size() == 0) {
        log_info("Diffusive ABM does not contain any agents. Population is initialized from SMM population.");
        int num_agents = 0;
        while (num_agents < total_pop) {
            auto position        = pos_rng();
            auto infection_state = state_rng(thread_local_rng(), current_pop);
            abm_pop.push_back(SingleWell<InfectionState>::Agent{position, InfectionState(infection_state)});
            current_pop[infection_state] = std::max(0., current_pop[infection_state] - 1);
            num_agents++;
        }
    }
    else {
        assert(int(total_pop) == int(abm_pop.size()) && "Population sizes of dabm and smm do not match.");
        for (auto& a : abm_pop) {
            auto infection_state         = state_rng(thread_local_rng(), current_pop);
            a.position                   = pos_rng();
            a.status                     = InfectionState(infection_state);
            current_pop[infection_state] = std::max(0., current_pop[infection_state] - 1);
        }
    }
}

} //namespace hybrid

} //namespace mio
