/*
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Anna Wendler
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
#ifndef IDE_INITIALFLOWS_H
#define IDE_INITIALFLOWS_H

#include "ide_secir/parameters.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/model.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"

namespace mio
{
namespace isecir
{
ScalarType compute_global_support_max(Model model, ScalarType dt);

void compute_previous_flows(Model& model, Eigen::Index idx_CurrentFlow, Eigen::Index idx_OutgoingFlow,
                            Eigen::Index time_series_index, ScalarType dt);

void compute_flows_with_mean(Model& model, Eigen::Index idx_CurrentFlow, Eigen::Index idx_OutgoingFlow, ScalarType dt,
                             Eigen::Index time_series_index);

// assume that we know flow from C to I from RKI data; compute the remaining flows based on this
void set_initial_flows(Model& model, ScalarType dt, ScalarType rki_cases_dummy, ScalarType rki_deaths_dummy);

void print_transitions(Model model);

} // namespace isecir
} // namespace mio

#endif // IDE_INITIALFLOWS_H