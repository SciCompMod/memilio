/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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

namespace mio
{
namespace abm
{

Simulation::Simulation(TimePoint t, World&& world)
    : m_world(std::move(world))
    , m_result(Eigen::Index(InfectionState::Count))
    , m_t(t)
    , m_dt(hours(1))
{
    initialize_locations(m_t);
}

void Simulation::initialize_locations(TimePoint t)
{
    for (auto& location : m_world.get_locations()) {
        location.initialize_subpopulations(t);
    }
}

void Simulation::advance(TimePoint tmax)
{
    //log initial system state
    initialize_locations(m_t);
    store_result_at(m_t);
    while (m_t < tmax) {
        evolve_world(tmax);
        store_result_at(m_t);
    }
}

void Simulation::evolve_world(TimePoint tmax)
{
    auto dt = std::min(m_dt, tmax - m_t);
    m_world.evolve(m_t, dt);
    m_t += m_dt;
}

void Simulation::store_result_at(TimePoint t)
{
    m_result.add_time_point(t.days());
    m_result.get_last_value().setZero();
    for (auto& location : m_world.get_locations()) {
        m_result.get_last_value() += location.get_subpopulations().get_last_value().cast<ScalarType>();
    }
}

} // namespace abm
} // namespace mio
