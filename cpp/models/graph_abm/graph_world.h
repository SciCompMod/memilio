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
#ifndef EPI_GRAPH_ABM_WORLD_H
#define EPI_GRAPH_ABM_WORLD_H

#include "abm/world.h"
#include "abm/location_type.h"

#include <vector>
#include <memory>

namespace mio
{
namespace graph_abm
{

class GraphWorld : public mio::abm::World
{
    //use abm::World constructors
    using mio::abm::World::World;
    using Base = mio::abm::World;

public:
    mio::abm::LocationId add_external_location(mio::abm::LocationType type, uint32_t num_cells)
    {
        mio::abm::LocationId id = {static_cast<uint32_t>(m_locations_external.size()), type};
        m_locations_external.emplace_back(std::make_unique<mio::abm::Location>(id, num_cells));
        return id;
    }
    mio::abm::Person& add_person(const mio::abm::LocationId id, mio::abm::AgeGroup age)
    {
        uint32_t person_id = static_cast<uint32_t>(m_persons_internal.size());
        m_persons_internal.push_back(
            std::make_shared<mio::abm::Person>(get_individualized_location(id), age, person_id));
        auto& person = *m_persons_internal.back();
        get_individualized_location(id).add_person(person);
        return person;
    }

private:
    std::vector<std::shared_ptr<mio::abm::Person>> m_persons_internal;
    std::vector<std::shared_ptr<mio::abm::Person>> m_persons_external;
    //std::vector<std::unique_ptr<mio::abm::Location>> m_locations_internal;
    std::vector<std::unique_ptr<mio::abm::Location>> m_locations_external;
    std::vector<std::pair<mio::abm::LocationId, mio::abm::LocationId>> m_external_location_mapping;
};

} // namespace graph_abm
} // namespace mio

#endif
