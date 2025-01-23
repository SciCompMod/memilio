/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Maximilian Betz
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
#ifndef PYMIO_SIMULATION_DAY_H
#define PYMIO_SIMULATION_DAY_H

#include "memilio/epidemiology/simulation_day.h"

#include "pybind11/pybind11.h"

namespace pymio
{

template <>
inline std::string pretty_name<mio::SimulationDay>()
{
    return "SimulationDay";
}

} // namespace pymio

#endif //PYMIO_SIMULATION_DAY_H
