/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Maximilian Betz
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
#ifndef PYMIO_OSEIR_H
#define PYMIO_OSEIR_H

//Includes from pymio
#include "pybind_util.h"
#include "utils/index.h"
#include "utils/custom_index_array.h"
#include "utils/parameter_set.h"
#include "compartments/simulation.h"
#include "compartments/flow_simulation.h"
#include "compartments/compartmentalmodel.h"
#include "epidemiology/age_group.h"
#include "epidemiology/populations.h"

//Includes from MEmilio
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "memilio/data/analyze_result.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;

namespace pymio
{
//specialization of pretty_name
template <>
inline std::string pretty_name<mio::oseir::InfectionState>()
{
    return "InfectionState";
}

} // namespace pymio

void bind_oseir(py::module_& m);

#endif //PYMIO_OSEIR_H
