/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: David Kermann, Sascha Korf
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
#ifndef EPI_ABM_TESTING_SET_H
#define EPI_ABM_TESTING_SET_H

namespace mio
{
namespace abm
{

class TestingScheme
{
public:
    TestingScheme() = default;

    TestingScheme(const std::vector<TestRule>& rules);

    add_rule(const TestRule& rule);

    run_scheme(const Person& person, const Time& time, const LocationType& location_type);

    void check_rules_for_activeness() const;

    bool get_active_status() const;
    {
        return m_active;
    }
    void set_active_status(const bool active)
    {
        m_active = active;
    }

private:
    std::vector<TestRule> m_test_rules = {};
    bool m_active                      = false;
}

class TestRule
{
public:
    TestRule() = default;
    TestRule(const std::vector<AgeGroup>& ageGroups, const std::vector<LocationType>& locationsGroups,
             const std::vector<InfectionState>& infectionGroups, TimeSpan interval, double probability);

    bool test_person(const Person& person) const; // Returns if the person's test.

    TimeSpan get_interval() const
    {
        return m_time_interval;
    }
    void set_interval(TimeSpan ts)
    {
        m_time_interval = t;
    }

    double get_probability() const
    {
        return m_probability;
    }
    void set_probability(double p)
    {
        m_probability = p;
    }

    bool active_status(TimePoint t) const
    {
        return;
    }

private:
    std::vector<AgeGroup> m_tested_ages                   = {};
    std::vector<LocationType> m_tested_locations          = {};
    std::vector<InfectionState> m_tested_infection_states = {};

    std::pair<TimeSpan, TimeSpan> m_time_interval; // Time Span where this test is active.
    double m_probability; // Probability that a given person actually tests themself.
    LocalTestingParameters m_testing_parameters; // Local testing parameters.
}

} // namespace abm
} // namespace mio

#endif