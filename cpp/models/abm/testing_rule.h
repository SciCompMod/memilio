/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: David Kerkmann
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
#ifndef EPI_ABM_TESTING_RULE_H
#define EPI_ABM_TESTING_RULE_H

#include "abm/person.h"
#include "abm/location.h"

namespace mio
{
namespace abm
{

class TestingRule
{
public:
    /**
     * Create a testing rule.
     * @param ages vector of age groups that are allowed for testing
     * @param location_types vector of location types that are allowed for testing
     * @param infection_states vector of infection states that are allowed for testing
     * An empty vector of ages/location types/infection states means that no condition on the corresponding property is set!
     */
    TestingRule(const std::vector<AgeGroup> ages, const std::vector<LocationType> location_types,
                const std::vector<InfectionState> infection_states);

    /**
     * Compares two testing rules for equality.
     * Compare references. Still possible to clone rules.
     */
    bool operator==(const TestingRule& other) const
    {
        return this == &other;
    }

    /**
     * add an age group to the set of age groups that are allowed for testing
     * @param age_group age group to be added
     */
    void add_age_group(const AgeGroup age_group);

    /**
     * remove an age group from the set of age groups that are allowed for testing
     * @param age_group age group to be removed
     */
    void remove_age_group(const AgeGroup age_group);

    /**
     * add a location type to the set of location types that are allowed for testing
     * @param location_type location type to be added
     */
    void add_location_type(const LocationType location_type);

    /**
     * remove a location tpye from the set of location tpyes that are allowed for testing
     * @param location_type location type to be removed
     */
    void remove_location_type(const LocationType location_type);

    /**
     * add an infection state to the set of infection states that are allowed for testing
     * @param infection_state infection state to be added
     */
    void add_infection_state(const InfectionState infection_state);

    /**
     * remove an infection state from the set of infection states that are allowed for testing
     * @param infection_state infection state to be removed
     */
    void remove_infection_state(const InfectionState infection_state);

    /**
     * check if a person and a location meet all the required properties to get tested
     * @param p person to be checked
     * @param l location to be checked
     */
    bool evaluate(const Person& p, const Location& l) const;

private:
    /**
     * check if a person has the required age to get tested
     * @param p person to be checked
     */
    bool has_requested_age(const Person& p) const;

    /**
     * check if a location is in the set of locations that are allowed for testing
     * @param l location to be checked
     */
    bool is_requested_location_type(const Location& l) const;

    /**
     * check if a person has the required infection state to get tested
     * @param p person to be checked
     */
    bool has_requested_infection_state(const Person& p) const;
    std::vector<AgeGroup> m_ages;
    std::vector<LocationType> m_location_types;
    std::vector<InfectionState> m_infection_states;
};

} // namespace abm
} // namespace mio

#endif /* testing_rule_h */
