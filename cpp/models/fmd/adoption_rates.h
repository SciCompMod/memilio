/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Kilian Volmer
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
#ifndef FMD_ADOPTION_RATES_H
#define FMD_ADOPTION_RATES_H

#include "memilio/config.h"
#include "memilio/geography/geolocation.h"
#include "memilio/geography/rtree.h"
#include "memilio/geography/distance.h"
#include "memilio/mobility/graph_simulation.h"
#include "memilio/mobility/metapopulation_mobility_asymmetric.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/graph_builder.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/timer/auto_timer.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/random_number_generator.h"
#include "smm/simulation.h"
#include "smm/parameters.h"
#include "fmd/infection_state.h"
#include <vector>

namespace mio
{
namespace fmd
{
using Status = mio::Index<InfectionState>;
using mio::regions::Region;

using AR = mio::AdoptionRate<ScalarType, Status, Region>;

std::vector<AR> generic_adoption_rates()
{
    auto home = Region(0);
    std::vector<AR> adoption_rates;
    adoption_rates.push_back({InfectionState::S,
                              InfectionState::E,
                              home,
                              0.2,
                              {{InfectionState::I, 0.8}, {InfectionState::INS, 0.1}, {InfectionState::ICS, 0.5}}});
    adoption_rates.push_back({InfectionState::E, InfectionState::I, home, 0.2, {}});
    adoption_rates.push_back({InfectionState::I, InfectionState::INS, home, 0.1, {}});
    adoption_rates.push_back({InfectionState::I, InfectionState::ICS, home, 0.1, {}});
    adoption_rates.push_back({InfectionState::ICS, InfectionState::D, home, 0.6, {}});
    adoption_rates.push_back({InfectionState::ICS, InfectionState::R, home, 0.4, {}});
    adoption_rates.push_back({InfectionState::INS, InfectionState::R, home, 0.5, {}});
    return adoption_rates;
}

} // namespace fmd
} // namespace mio

#endif // FMD_ADOPTION_RATES_H