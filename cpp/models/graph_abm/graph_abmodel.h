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

#ifndef MIO_ABM_GRAPH_ABMODEL_H
#define MIO_ABM_GRAPH_ABMODEL_H

#include "abm/location_type.h"
#include "abm/model.h"
#include "abm/person_id.h"
#include "abm/time.h"
#include "abm/location_id.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/mioomp.h"
#include "abm/mobility_rules.h"
#include "abm/mobility_rules.h"
#include <cstddef>
#include <cstdint>
#include <list>
#include <vector>

namespace mio
{
using namespace abm;
class GraphABModel : public abm::Model
{
    using Base = Model;

public:
    GraphABModel(size_t num_agegroups, int id,
                 std::vector<Base::MobilityRuleType> mobility_rules =
                     std::vector<Base::MobilityRuleType>{&get_buried, &return_home_when_recovered, &go_to_hospital,
                                                         &go_to_icu, &go_to_school, &go_to_work, &go_to_shop,
                                                         &go_to_event, &go_to_quarantine})
        : Base(num_agegroups, id)
    {
        Base::m_mobility_rules = mobility_rules;
    }

    /**
     * @brief Get person buffer. 
     */
    std::vector<size_t>& get_person_buffer()
    {
        return m_person_buffer;
    }

    /** 
    * @brief Removes person from the model.
    * @param[in] pos Index of person in m_persons vector.
    */
    void remove_person(size_t pos)
    {
        Base::m_persons.erase(Base::m_persons.begin() + pos);
        Base::m_activeness_statuses.erase(Base::m_activeness_statuses.begin() + pos);
        if (Base::m_person_ids_equal_index) {
            Base::m_person_ids_equal_index = false;
        }
    }

    /** 
     * @brief Evolve the Graph Model one time step.
     * @param[in] t Current time.
     * @param[in] dt Length of the time step.
     */
    void evolve(TimePoint t, TimeSpan dt)
    {
        Base::begin_step(t, dt);
        log_info("Graph ABM Model interaction.");
        Base::interaction(t, dt);
        log_info("Graph ABM Model mobility.");
        perform_mobility(t, dt);
    }

private:
    void perform_mobility(TimePoint t, TimeSpan dt)
    {
        const uint32_t num_persons = static_cast<uint32_t>(Base::m_persons.size());
        for (uint32_t person_index = 0; person_index < num_persons; ++person_index) {
            if (Base::m_activeness_statuses[person_index]) {
                Person& person    = Base::m_persons[person_index];
                auto personal_rng = PersonalRandomNumberGenerator(Base::m_rng, person);

                auto try_mobility_rule = [&](auto rule) -> bool {
                    //run mobility rule and check if change of location can actually happen
                    auto target_type = rule(personal_rng, person, t, dt, parameters);
                    if (person.get_assigned_location_model_id(target_type) == Base::m_id) {
                        const Location& target_location = Base::get_location(Base::find_location(target_type, person));
                        const LocationId current_location = person.get_location();
                        // the Person cannot move if they do not wear mask as required at targeted location
                        if (target_location.is_mask_required() &&
                            !person.is_compliant(personal_rng, InterventionType::Mask)) {
                            return false;
                        }
                        // the Person cannot move if the capacity of targeted Location is reached
                        if (target_location.get_id() == current_location ||
                            get_number_persons(target_location.get_id()) >= target_location.get_capacity().persons) {
                            return false;
                        }
                        // the Person cannot move if the performed TestingStrategy is positive
                        if (!m_testing_strategy.run_strategy(personal_rng, person, target_location, t)) {
                            return false;
                        }
                        // update worn mask to target location's requirements
                        if (target_location.is_mask_required()) {
                            // if the current MaskProtection level is lower than required, the Person changes mask
                            if (parameters.get<MaskProtection>()[person.get_mask().get_type()] <
                                parameters.get<MaskProtection>()[target_location.get_required_mask()]) {
                                person.set_mask(target_location.get_required_mask(), t);
                            }
                        }
                        else {
                            person.set_mask(MaskType::None, t);
                        }
                        Base::change_location(person, target_location.get_id());
                        return true;
                    }
                    else { //person moves to other world
                        Base::m_activeness_statuses[person_index] = false;
                        person.set_location(target_type, abm::LocationId::invalid_id(),
                                            std::numeric_limits<int>::max());
                        m_person_buffer.push_back(person_index);
                        m_are_exposure_caches_valid       = false;
                        m_is_local_population_cache_valid = false;
                        return true;
                    }
                };

                for (auto rule : Base::m_mobility_rules) {
                    bool applied = try_mobility_rule(rule);
                    //only use one mobility rule per person
                    if (applied) {
                        break;
                    }
                }
            }
        }

        // check if a person makes a trip
        bool weekend     = t.is_weekend();
        size_t num_trips = Base::m_trip_list.num_trips(weekend);

        for (; Base::m_trip_list.get_current_index() < num_trips &&
               Base::m_trip_list.get_next_trip_time(weekend).seconds() < (t + dt).time_since_midnight().seconds();
             Base::m_trip_list.increase_index()) {
            auto& trip        = Base::m_trip_list.get_next_trip(weekend);
            auto& person      = get_person(trip.person_id);
            auto person_index = Base::get_person_index(trip.person_id);
            auto personal_rng = PersonalRandomNumberGenerator(m_rng, person);
            // skip the trip if the person is in quarantine or is dead
            if (person.is_in_quarantine(t, parameters) || person.get_infection_state(t) == InfectionState::Dead) {
                continue;
            }
            if (trip.destination_model_id == Base::m_id) {
                auto& target_location = get_location(trip.destination);
                // skip the trip if the Person wears mask as required at targeted location
                if (target_location.is_mask_required() && !person.is_compliant(personal_rng, InterventionType::Mask)) {
                    continue;
                }
                // skip the trip if the performed TestingStrategy is positive
                if (!Base::m_testing_strategy.run_strategy(personal_rng, person, target_location, t)) {
                    continue;
                }
                // all requirements are met, move to target location
                change_location(person, target_location.get_id(), trip.trip_mode);
                // update worn mask to target location's requirements
                if (target_location.is_mask_required()) {
                    // if the current MaskProtection level is lower than required, the Person changes mask
                    if (parameters.get<MaskProtection>()[person.get_mask().get_type()] <
                        parameters.get<MaskProtection>()[target_location.get_required_mask()]) {
                        person.set_mask(target_location.get_required_mask(), t);
                    }
                }
                else {
                    person.set_mask(MaskType::None, t);
                }
            }
            else { //person moves to other world
                Base::m_activeness_statuses[person_index] = false;
                person.set_location(trip.destination_type, abm::LocationId::invalid_id(),
                                    std::numeric_limits<int>::max());
                m_person_buffer.push_back(person_index);
                m_are_exposure_caches_valid       = false;
                m_is_local_population_cache_valid = false;
            }
        }
        if (((t).days() < std::floor((t + dt).days()))) {
            Base::m_trip_list.reset_index();
        }
    }

    std::vector<size_t> m_person_buffer; ///< List with indices of persons that are subject to move to another node.
};
} // namespace mio

#endif //MIO_ABM_GRAPH_ABMODEL_H
