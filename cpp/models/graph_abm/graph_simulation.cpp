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
    Base::initialize_locations(t);
    Base::store_result_at(t);
}

} // namespace graph_abm
} // namespace mio
