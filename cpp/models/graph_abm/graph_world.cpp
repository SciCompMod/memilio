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
#include "models/abm/world.h"
#include "graph_abm/graph_world.h"

namespace mio
{
namespace graph_abm
{

mio::abm::Location& GraphWorld::find_location(mio::abm::LocationType type, const mio::abm::Person& person)
{
    auto index    = person.get_assigned_location_index(type);
    auto world_id = person.get_assigned_location_world_id(type);
    return get_individualized_location({index, type, world_id});
}

mio::abm::Location& GraphWorld::get_individualized_location(mio::abm::LocationId id)
{
    if (id.world_id == Base::m_world_id) {
        return *Base::m_locations[id.index];
    }
    else {
        mio::abm::Location loc = mio::abm::Location(id.type, id.index, id.world_id);
        auto iter              = std::find_if(m_locations_external.begin(), m_locations_external.end(),
                                              [loc](const std::unique_ptr<mio::abm::Location>& loc_ext) {
                                     return ((loc.get_type() == loc_ext->get_type()) &&
                                             (loc.get_index() == loc_ext->get_index()) &&
                                             (loc.get_world_id() == loc_ext->get_world_id()));
                                 });
        if (iter != m_locations_external.end()) {
            return *(*iter);
        }
        else {
            m_locations_external.emplace_back(std::make_unique<mio::abm::Location>(loc));
            return *m_locations_external.back();
        }
    }
}

void GraphWorld::add_existing_person(std::unique_ptr<mio::abm::Person>&& person)
{
    Base::m_persons.push_back(std::move(person));
}

std::unique_ptr<mio::abm::Person>& GraphWorld::get_person(uint32_t person_id, uint32_t person_world_id)
{
    return *std::find_if(Base::m_persons.begin(), Base::m_persons.end(),
                         [person_id, person_world_id](const std::unique_ptr<mio::abm::Person>& person) {
                             return (person_id == person->get_person_id() && person_world_id == person->get_world_id());
                         });
}

void GraphWorld::migration(mio::abm::TimePoint t, mio::abm::TimeSpan dt)
{
    size_t person_iter = 0;
    while (person_iter < Base::m_persons.size()) {
        for (auto rule : Base::m_migration_rules) {
            //check if transition rule can be applied
            auto target_type      = rule.first(*Base::m_persons[person_iter], t, dt, m_migration_parameters);
            auto& target_location = find_location(target_type, *Base::m_persons[person_iter]);
            auto current_location = Base::m_persons[person_iter]->get_location();
            if (m_testing_strategy.run_strategy(*Base::m_persons[person_iter], target_location, t)) {
                if (target_location != current_location &&
                    target_location.get_number_persons() < target_location.get_capacity().persons) {
                    bool wears_mask = Base::m_persons[person_iter]->apply_mask_intervention(target_location);
                    if (wears_mask) {
                        if (target_location.get_world_id() != Base::m_world_id) {
                            //person changes world
                            Base::m_persons[person_iter]->migrate_to_other_world(target_location, false);
                            m_persons_to_migrate.push_back(std::move(Base::m_persons[person_iter]));
                            Base::m_persons.erase(Base::m_persons.begin() + person_iter);
                            goto CONTINUE_PERSON;
                        }
                        else {
                            Base::m_persons[person_iter]->migrate_to(target_location);
                        }
                    }
                }
            }
        }
        ++person_iter;
    CONTINUE_PERSON:
        continue;
    }
    // check if a person makes a trip
    size_t num_trips = m_trip_list.num_trips();
    if (num_trips != 0) {
        while (m_trip_list.get_current_index() < num_trips && m_trip_list.get_next_trip_time() < t + dt) {
            auto& trip            = m_trip_list.get_next_trip();
            auto& person          = get_person(trip.person_id, trip.person_world_id);
            auto current_location = person->get_location();
            if (!person->is_in_quarantine() && current_location == get_individualized_location(trip.migration_origin)) {
                auto& target_location = get_individualized_location(trip.migration_destination);
                person->apply_mask_intervention(target_location);
                if (m_testing_strategy.run_strategy(*person, target_location, t)) {
                    if (target_location.get_world_id() != Base::m_world_id) {
                        //person changes world
                        person->migrate_to_other_world(target_location, false);
                        m_persons_to_migrate.push_back(std::move(person));
                        Base::m_persons.erase(Base::m_persons.begin() + trip.person_id);
                    }
                    else {
                        person->migrate_to(target_location);
                    }
                }
            }
            m_trip_list.increase_index();
        }
    }
}

std::vector<std::unique_ptr<mio::abm::Person>>& GraphWorld::get_persons_to_migrate()
{
    return m_persons_to_migrate;
}

void GraphWorld::evolve(mio::abm::TimePoint t, mio::abm::TimeSpan dt)
{
    Base::begin_step(t, dt);
    Base::interaction(t, dt);
    Base::get_testing_strategy().update_activity_status(t);
    migration(t, dt);
    Base::end_step(t, dt);
}

} // namespace graph_abm
} // namespace mio
