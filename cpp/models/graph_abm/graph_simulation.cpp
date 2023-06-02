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
#include "graph_abm/graph_simulation.h"

namespace mio
{
namespace graph_abm
{
GraphSimulation::GraphSimulation(mio::abm::TimePoint t, GraphWorld&& graph_world)
    : Base(t)
    , m_graph_world(std::move(graph_world))
{
    Base::m_dt     = mio::abm::hours(1);
    Base::m_result = mio::TimeSeries<ScalarType>(Eigen::Index(mio::abm::InfectionState::Count));
    initialize_locations(t);
    store_result_at(t);
}

void GraphSimulation::initialize_locations(mio::abm::TimePoint t)
{
    for (auto& location : m_graph_world.get_locations()) {
        location.initialize_subpopulations(t);
    }
}

void GraphSimulation::store_result_at(mio::abm::TimePoint t)
{
    Base::m_result.add_time_point(t.days());
    Base::m_result.get_last_value().setZero();
    for (auto& location : m_graph_world.get_locations()) {
        Base::m_result.get_last_value() += location.get_subpopulations().get_last_value().cast<ScalarType>();
    }
}

void GraphSimulation::advance(mio::abm::TimePoint tmax)
{
    auto t = Base::m_t;
    while (t < tmax) {
        auto dt = std::min(Base::m_dt, tmax - t);
        m_graph_world.evolve(t, dt);
        t += Base::m_dt;
        store_result_at(t);
    }
    Base::m_t = tmax;
}

} // namespace graph_abm
} // namespace mio
