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
#include "abm/simulation.h"

namespace mio
{

AbmSimulation::AbmSimulation(TimePoint t, World&& world)
    : m_world(std::move(world))
    , m_result(Eigen::Index(InfectionState::Count))
    , m_t(t)
    , m_dt(hours(1))
{
    store_result_at(t);
}

void AbmSimulation::advance(TimePoint tmax)
{
    auto t = m_t;
    while (t < tmax) {
        auto dt = std::min(m_dt, tmax - t);
        m_world.evolve(t, dt);
        t += m_dt;
        store_result_at(t);
    }
}

void AbmSimulation::store_result_at(TimePoint t)
{
    m_result.add_time_point(t.days());
    m_result.get_last_value().setZero();
    for (auto&& locations : m_world.get_locations()) {
        for (auto& location : locations){
            m_result.get_last_value() += location.get_subpopulations().cast<double>();
        }
    }
}

} // namespace mio