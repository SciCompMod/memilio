/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Julia Bicker
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

#ifndef MIO_ABM_MODEL_WRAPPER_H
#define MIO_ABM_MODEL_WRAPPER_H

#include "abm/model.h"
#include "abm/time.h"
#include "abm/location_id.h"
#include "memilio/utils/mioomp.h"
#include "abm/mobility_rules.h"
#include <list>

namespace mio
{
using namespace abm;
class ModelWrapper : public abm::Model
{
    using Model::Model;
    using Base = Model;

public:
    /**
     * @brief Get a reference to a Person from this Model.
     * @param[in] id A person's PersonId. 
     * First 32 bit are the Person's index and second 32 bit the Persons's home model id. 
     * @return Position of Person in m_persons vector.
     * @{
     */
    size_t get_person_pos(PersonId id)
    {
        auto it = std::find(Base::m_persons.begin(), Base::m_persons.end(), [id](auto& person) {
            person.get_id() == id;
        });
        if (it == Base::m_persons.end()) {
            log_error("Given PersonId is not in this Model.");
            return std::numeric_limits<size_t>::max();
        }
        else {
            return std::distance(Base::m_persons.begin(), it);
        }
    }

private:
    void perform_mobility(TimePoint t, TimeSpan dt)
    {
        const uint32_t num_persons = static_cast<uint32_t>(Base::m_persons.size());
    PRAGMA_OMP(parallel for)
    for (uint32_t person_id = 0; person_id < num_persons; ++person_id) {
        if (Base::m_activeness_statuses[person_id]) {
            Person& person    = Base::m_persons[person_id];
            auto personal_rng = PersonalRandomNumberGenerator(Base::m_rng, person);

            auto try_mobility_rule = [&](auto rule) -> bool {
                //run mobility rule and check if change of location can actually happen
                auto target_type = rule(personal_rng, person, t, dt, parameters);
                if (person.get_assigned_location_model_id(target_type) == Base::m_id) {
                    const Location& target_location   = Base::get_location(Base::find_location(target_type, person_id));
                    const LocationId current_location = person.get_location();
                    if (Base::m_testing_strategy.run_strategy(personal_rng, person, target_location, t)) {
                        if (target_location.get_id() != current_location &&
                            Base::get_number_persons(target_location.get_id()) <
                                target_location.get_capacity().persons) {
                            bool wears_mask = person.apply_mask_intervention(personal_rng, target_location);
                            if (wears_mask) {
                                Base::change_location(person_id, target_location.get_id());
                            }
                            return true;
                        }
                    }
                }
                else { //person moves to other world
                    Base::m_activeness_statuses[person_id] = false;
                    person.set_location(target_type, abm::LocationId::invalid_id(), std::numeric_limits<int>::max());
                    person_buffer.push_back(person_id);
                    return true;
                }
                return false;
            };

            //run mobility rules one after the other if the corresponding location type exists
            //shortcutting of bool operators ensures the rules stop after the first rule is applied
            if (m_use_mobility_rules) {
                (Base::has_locations({LocationType::Cemetery}) && try_mobility_rule(&get_buried)) ||
                    (Base::has_locations({LocationType::Home}) && try_mobility_rule(&return_home_when_recovered)) ||
                    (Base::has_locations({LocationType::Hospital}) && try_mobility_rule(&go_to_hospital)) ||
                    (Base::has_locations({LocationType::ICU}) && try_mobility_rule(&go_to_icu)) ||
                    (Base::has_locations({LocationType::School, LocationType::Home}) &&
                     try_mobility_rule(&go_to_school)) ||
                    (Base::has_locations({LocationType::Work, LocationType::Home}) && try_mobility_rule(&go_to_work)) ||
                    (Base::has_locations({LocationType::BasicsShop, LocationType::Home}) &&
                     try_mobility_rule(&go_to_shop)) ||
                    (Base::has_locations({LocationType::SocialEvent, LocationType::Home}) &&
                     try_mobility_rule(&go_to_event)) ||
                    (Base::has_locations({LocationType::Home}) && try_mobility_rule(&go_to_quarantine));
            }
            else {
                //no daily routine mobility, just infection related
                (Base::has_locations({LocationType::Cemetery}) && try_mobility_rule(&get_buried)) ||
                    (Base::has_locations({LocationType::Home}) && try_mobility_rule(&return_home_when_recovered)) ||
                    (Base::has_locations({LocationType::Hospital}) && try_mobility_rule(&go_to_hospital)) ||
                    (Base::has_locations({LocationType::ICU}) && try_mobility_rule(&go_to_icu)) ||
                    (Base::has_locations({LocationType::Home}) && try_mobility_rule(&go_to_quarantine));
            }
        }
    }

    // check if a person makes a trip
    bool weekend     = t.is_weekend();
    size_t num_trips = m_trip_list.num_trips(weekend);

    if (num_trips != 0) {
        while (m_trip_list.get_current_index() < num_trips &&
               m_trip_list.get_next_trip_time(weekend).seconds() < (t + dt).time_since_midnight().seconds()) {
            auto& trip        = m_trip_list.get_next_trip(weekend);
            auto person_pos   = get_person_pos(trip.person_id);
            auto& person      = Base::get_person(person_pos);
            auto personal_rng = PersonalRandomNumberGenerator(m_rng, person);
            if (!person.is_in_quarantine(t, parameters) && person.get_infection_state(t) != InfectionState::Dead) {
                if (trip.destination_model_id == Base::m_id) {
                    auto& target_location = get_location(trip.destination);
                    if (m_testing_strategy.run_strategy(personal_rng, person, target_location, t)) {
                        person.apply_mask_intervention(personal_rng, target_location);
                        change_location(person.get_id(), target_location.get_id(), trip.trip_mode);
                    }
                }
                else {
                    //person moves to other world
                    Base::m_activeness_statuses[person_pos] = false;
                    person.set_location(trip.destination_type, abm::LocationId::invalid_id(),
                                        std::numeric_limits<int>::max());
                    person_buffer.push_back(person_pos);
                }
            }
            m_trip_list.increase_index();
        }
    }
    if (((t).days() < std::floor((t + dt).days()))) {
        m_trip_list.reset_index();
    }
    }

    std::list<size_t> person_buffer; ///< List with indices of persons that are deactivated.
};
} // namespace mio

#endif //MIO_ABM_MODEL_WRAPPER_H