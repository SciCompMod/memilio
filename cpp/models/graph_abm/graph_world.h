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
#include <map>

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
    // mio::abm::LocationId add_external_location(mio::abm::LocationType type, uint32_t num_cells = 1)
    // {
    //     mio::abm::LocationId id = {static_cast<uint32_t>(m_locations_external.size() + Base::m_locations.size()), type};
    //     m_locations_external.emplace_back(std::make_unique<mio::abm::Location>(id, num_cells));
    //     return id;
    // }

    // mio::abm::Person& add_person(const mio::abm::LocationId id, mio::abm::AgeGroup age)
    // {
    //     uint32_t person_id = static_cast<uint32_t>(Base::m_persons.size());
    //     Base::m_persons.push_back(
    //         std::make_unique<mio::abm::Person>(get_individualized_location(id), age, person_id, m_world_id));
    //     auto& person = *m_persons_internal.back();
    //     get_individualized_location(id).add_person(person);
    //     return person;
    // }

    // void insert_external_location_to_map(int node_to, mio::abm::LocationId id_this_world,
    //                                      mio::abm::LocationId id_other_world)
    // {
    //     if (m_external_location_mapping.find(node_to) != m_external_location_mapping.end()) {
    //         auto& mapping = m_external_location_mapping[node_to];
    //         mapping.push_back(std::make_pair(id_this_world, id_other_world));
    //     }
    //     else {
    //         std::vector<std::pair<mio::abm::LocationId, mio::abm::LocationId>> v = {
    //             std::make_pair(id_this_world, id_other_world)};
    //         m_external_location_mapping.insert(
    //             std::pair<int, std::vector<std::pair<mio::abm::LocationId, mio::abm::LocationId>>>(node_to, v));
    //     }
    // }
    // auto get_persons() const -> Range<std::pair<ConstPersonIterator, ConstPersonIterator>>
    // {
    //     return std::make_pair(ConstPersonIterator(m_persons_internal.begin()),
    //                           ConstPersonIterator(m_persons_internal.end()));
    // }

    // /**
    //  * Prepare the World for the next simulation step.
    //  * @param[in] t Current time.
    //  * @param[in] dt Length of the time step.
    //  */
    // void begin_step(mio::abm::TimePoint t, mio::abm::TimeSpan dt);

    /** @brief Person%s move in the World according to rules.
     * @param[in] t The current TimePoint.
     * @param[in] dt The length of the time step of the Simulation.
     */
    void migration(mio::abm::TimePoint t, mio::abm::TimeSpan dt);

    /**
     * Evolve the world one time step.
     * @param[in] t Current time.
     * @param[in] dt Length of the time step.
     */
    void evolve(mio::abm::TimePoint t, mio::abm::TimeSpan dt);

        /**
     * @brief Find an assigned Location of a Person.
     * @param[in] type The #LocationType that specifies the assigned Location.
     * @param[in] person The Person.
     * @return Pointer to the assigned Location. If the location is in another world, the function will return the location in m_locations_external.
     */
    mio::abm::Location& find_location(mio::abm::LocationType type, const mio::abm::Person& person);

private:
    // /**
    //  * @brief Person%s interact at their Location and may become infected.
    //  * @param[in] t The current TimePoint.
    //  * @param[in] dt The length of the time step of the Simulation.
    //  */
    // void interaction(mio::abm::TimePoint t, mio::abm::TimeSpan dt);

    //persons that will migrate to other worlds
    std::vector<std::unique_ptr<mio::abm::Person>> m_persons_to_migrate;
    std::vector<std::unique_ptr<mio::abm::Location>> m_locations_external;
};

} // namespace graph_abm
} // namespace mio

#endif
