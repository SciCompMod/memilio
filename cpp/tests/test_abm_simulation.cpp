/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth, David Kerkmann, Sascha Korf, Martin J. Kuehn, Khoa Nguyen
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
#include "test_abm.h"

TEST(TestSimulation, advance_random)
{
    auto world     = mio::abm::World();
    auto location1 = world.add_location(mio::abm::LocationType::School);
    auto location2 = world.add_location(mio::abm::LocationType::School);
    auto& p1       = world.add_person(location1, mio::abm::InfectionState::Carrier, mio::abm::AgeGroup::Age5to14);
    auto& p2       = world.add_person(location1, mio::abm::InfectionState::Susceptible, mio::abm::AgeGroup::Age5to14);
    auto& p3       = world.add_person(location2, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age5to14);
    auto& p4       = world.add_person(location2, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age5to14);
    p1.set_assigned_location(location1);
    p2.set_assigned_location(location1);
    p3.set_assigned_location(location2);
    p4.set_assigned_location(location2);

    auto sim = mio::abm::Simulation(mio::abm::TimePoint(0), std::move(world));

    sim.advance(mio::abm::TimePoint(0) + mio::abm::hours(50));
    ASSERT_EQ(sim.get_result().get_num_time_points(), 51);
    ASSERT_THAT(sim.get_result().get_times(), ElementsAreLinspace(0.0, 50.0 / 24.0, 51));
    for (auto&& v : sim.get_result()) {
        ASSERT_EQ(v.sum(), 4);
    }
}

TEST(TestDiscreteDistribution, generate)
{
    using namespace mio;
    auto distribution = mio::DiscreteDistribution<size_t>();

    std::vector<double> weights;
    for (size_t i = 0; i < 50; i++) {
        weights = {};
        ASSERT_EQ(distribution(weights), 0);

        weights = {0.5};
        ASSERT_EQ(distribution(weights), 0);

        weights = {0.5, 1.3, 0.1, 0.4, 0.3};
        auto d  = distribution(weights);
        ASSERT_GE(d, 0);
        ASSERT_LE(d, 4);
    }
}

TEST(TestSimulation, advance_subpopulation)
{
    auto world       = mio::abm::World();
    auto location_id = world.add_location(mio::abm::LocationType::School);
    auto& school     = world.get_individualized_location(location_id);
    auto person1 = mio::abm::Person(location_id, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age5to14, {});
    school.add_person(person1);
    auto person2 = mio::abm::Person(location_id, mio::abm::InfectionState::Infected, mio::abm::AgeGroup::Age15to34, {});
    school.add_person(person2);
    auto person3 = mio::abm::Person(location_id, mio::abm::InfectionState::Exposed, mio::abm::AgeGroup::Age35to59, {});
    school.add_person(person3);

    auto sim = mio::abm::Simulation(mio::abm::TimePoint(0), std::move(world));
    sim.advance(mio::abm::TimePoint(0) + mio::abm::hours(50));

    for (size_t i = 0; i < 50; i++) {
        auto v = school.get_population().get_value(i);
        // Check whether the number of persons in infected state at the location is consistent
        ASSERT_LE(v[size_t(mio::abm::InfectionState::Infected)], 3);
        // Check the time evolution is correct
        ASSERT_EQ(school.get_population().get_time(i), ScalarType(i) / 24);
    }
}
