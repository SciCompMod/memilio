/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Khoa Nguyen
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
#include "abm/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/mioomp.h"
#include <random>

namespace mio
{
namespace abm
{

Simulation::Simulation(TimePoint t, Model&& model)
    : m_model(std::move(model))
    , m_t(t)
    , m_dt(hours(1))
{
}

void Simulation::evolve_model(TimePoint tmax)
{
    auto dt = std::min(m_dt, tmax - m_t);
    m_model.evolve(m_t, dt);
    m_t += m_dt;
}

} // namespace abm
} // namespace mio
