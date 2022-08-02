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
    TestingRule() = default;
    TestingRule(const std::vector<AgeGroup> ages, const std::vector<LocationType> location_type,
                const std::vector<InfectionState> infection_states);

    bool operator==(const TestingRule& other) const
    {
        return this == &other; // compare pointers. Still possible to clone Rules.
    }
    // nothing means everything is accepted!
    void add_age_group(const AgeGroup ageGroup);
    void remove_age_group(const AgeGroup ageGroup);

    void add_location_type(const LocationType locationType);
    void remove_location_type(const LocationType locationType);

    void add_infection_state(const InfectionState infection_state);
    void remove_infection_state(const InfectionState infection_state);

    bool evaluate(const Person& p, const Location& l) const;

private:
    bool has_requested_age(const Person& p) const;
    bool is_requested_location_type(const Location& l) const;
    bool has_requested_infection_state(const Person& p) const;
    std::vector<AgeGroup> m_ages                   = {};
    std::vector<LocationType> m_location_types     = {};
    std::vector<InfectionState> m_infection_states = {};
};

} // namespace abm
} // namespace mio

#endif /* testing_rule_h */
