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

#include "models/abm/world.h"
#include "models/abm/location_type.h"
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
    /**
     * @brief Get reference to all Person%s that will migrate to another world.
     * @return A reference to the vector of all Person%s that will migrate.
     */
    std::vector<std::unique_ptr<mio::abm::Person>>& get_persons_to_migrate();

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

    /**
     * @brief Add an existing Person to a world.
     * @param[in] person Unique pointer to the person that should beadded to world.
    */
    void add_existing_person(std::unique_ptr<mio::abm::Person>&& person);

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
