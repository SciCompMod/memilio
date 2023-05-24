/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: Julia Bicker, Daniel Abele, Martin J. Kuehn
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
#include "abm/world.h"
#include "graph_abm/graph_world.h"

namespace mio
{
namespace graph_abm
{

void GraphWorld::interaction(mio::abm::TimePoint t, mio::abm::TimeSpan dt)
{
    for (auto&& person : m_persons_internal) {
        person->interact(t, dt, m_infection_parameters);
    }
    for (auto&& person : m_persons_external) {
        person->interact(t, dt, m_infection_parameters);
    }
}

} // namespace graph_abm
} // namespace mio
