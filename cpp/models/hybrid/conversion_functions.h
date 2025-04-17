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

#ifndef MIO_HYBRID_CONVERSION_FUNCTIONS_H
#define MIO_HYBRID_CONVERSION_FUNCTIONS_H

#include "hybrid/temporal_hybrid_model.h"
#include "hybrid/infection_state.h"
#include "d_abm/simulation.h"
#include "d_abm/single_well.h"
#include "smm/simulation.h"

namespace mio
{
namespace hybrid
{

template <>
void convert_model(const dabm::Simulation<SingleWell<InfectionState>>& current_model,
                   smm::Simulation<1, InfectionState>& target_model);

template <>
void convert_model(const smm::Simulation<1, InfectionState>& current_model,
                   dabm::Simulation<SingleWell<InfectionState>>& target_model);

} //namespace hybrid

} //namespace mio

#endif //MIO_HYBRID_CONVERSION_FUNCTIONS_H
