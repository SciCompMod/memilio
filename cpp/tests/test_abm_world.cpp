/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
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

TEST(TestWorld, init)
{
    auto world = mio::abm::World(6);
    for (uint32_t i = 0; i < (uint32_t)mio::abm::LocationType::Count; i++) {
        ASSERT_THAT(world.get_locations()[i], testing::ElementsAre());
    }
    ASSERT_THAT(world.get_persons(), testing::ElementsAre());
}

TEST(TestWorld, copyWorld)
{
    auto simulation_params                              = mio::abm::SimulationParameters(6);
    simulation_params.get<mio::abm::IncubationPeriod>() = 4.;
    auto world                                          = mio::abm::World(simulation_params);
    world.add_location(mio::abm::LocationType::School);
    world.add_location(mio::abm::LocationType::School);
    world.add_location(mio::abm::LocationType::Work);
    world.add_location(mio::abm::LocationType::Home);

    auto copied_world = mio::abm::World(world);
    auto copied_simulation_params =
        copied_world.get_simulation_parameters()
            .get<mio::abm::IncubationPeriod>()[{mio::AgeGroup(0), mio::abm::VaccinationState::Count}];

    ASSERT_EQ(copied_world.get_locations().size(), (uint32_t)mio::abm::LocationType::Count);
    ASSERT_EQ(copied_world.get_locations()[(uint32_t)mio::abm::LocationType::School].size(), 2);
    ASSERT_EQ(copied_world.get_locations()[(uint32_t)mio::abm::LocationType::Work].size(), 1);
    ASSERT_EQ(copied_world.get_locations()[(uint32_t)mio::abm::LocationType::Home].size(), 1);
    ASSERT_EQ(copied_simulation_params.value(), 4.);
}

TEST(TestWorld, addLocation)
{
    auto world      = mio::abm::World(6);
    auto school_id1 = world.add_location(mio::abm::LocationType::School);
    auto school_id2 = world.add_location(mio::abm::LocationType::School);
    auto work_id    = world.add_location(mio::abm::LocationType::Work);
    auto home_id    = world.add_location(mio::abm::LocationType::Home);

    ASSERT_EQ((int)school_id1.index, 0);
    ASSERT_EQ((int)school_id2.index, 1);

    auto& school1 = world.get_individualized_location(school_id1);
    auto& school2 = world.get_individualized_location(school_id2);
    auto& work    = world.get_individualized_location(work_id);
    auto& home    = world.get_individualized_location(home_id);

    ASSERT_EQ(world.get_locations().size(), (uint32_t)mio::abm::LocationType::Count);
    ASSERT_EQ(world.get_locations()[(uint32_t)mio::abm::LocationType::School].size(), 2);

    ASSERT_EQ(&world.get_locations()[(uint32_t)mio::abm::LocationType::School][0], &school1);
    ASSERT_EQ(&world.get_locations()[(uint32_t)mio::abm::LocationType::School][1], &school2);
    ASSERT_EQ(&world.get_locations()[(uint32_t)mio::abm::LocationType::Work][0], &work);
    ASSERT_EQ(&world.get_locations()[(uint32_t)mio::abm::LocationType::Home][0], &home);
}

TEST(TestWorld, addPerson)
{
    auto world    = mio::abm::World(6);
    auto location = world.add_location(mio::abm::LocationType::School);

    auto& p1 = world.add_person(location, mio::abm::InfectionState::Recovered_Carrier);
    auto& p2 = world.add_person(location, mio::abm::InfectionState::Exposed);

    ASSERT_EQ(world.get_persons().size(), 2);
    ASSERT_EQ(&world.get_persons()[0], &p1);
    ASSERT_EQ(&world.get_persons()[1], &p2);
}

TEST(TestWorld, getSubpopulationCombined)
{
    auto world   = mio::abm::World(6);
    auto school1 = world.add_location(mio::abm::LocationType::School);
    auto school2 = world.add_location(mio::abm::LocationType::School);
    auto school3 = world.add_location(mio::abm::LocationType::School);
    world.add_person(school1, mio::abm::InfectionState::Carrier);
    world.add_person(school1, mio::abm::InfectionState::Susceptible);
    world.add_person(school2, mio::abm::InfectionState::Susceptible);
    world.add_person(school2, mio::abm::InfectionState::Susceptible);
    world.add_person(school3, mio::abm::InfectionState::Carrier);

    ASSERT_EQ(world.get_subpopulation_combined(mio::abm::InfectionState::Susceptible, mio::abm::LocationType::School),
              3);
    ASSERT_EQ(world.get_subpopulation_combined(mio::abm::InfectionState::Carrier, mio::abm::LocationType::School), 2);
}

TEST(TestWorld, findLocation)
{
    auto world     = mio::abm::World(6);
    auto home_id   = world.add_location(mio::abm::LocationType::Home);
    auto school_id = world.add_location(mio::abm::LocationType::School);
    auto work_id   = world.add_location(mio::abm::LocationType::Work);
    auto person    = mio::abm::Person(home_id, mio::abm::InfectionState::Recovered_Carrier, mio::AgeGroup(4),
                                      world.get_global_infection_parameters());
    auto& home     = world.get_individualized_location(home_id);
    auto& school   = world.get_individualized_location(school_id);
    auto& work     = world.get_individualized_location(work_id);
    person.set_assigned_location(home);
    person.set_assigned_location(school);
    person.set_assigned_location({0, mio::abm::LocationType::Work});

    ASSERT_EQ(world.find_location(mio::abm::LocationType::Work, person), &work);
    ASSERT_EQ(world.find_location(mio::abm::LocationType::School, person), &school);
    ASSERT_EQ(world.find_location(mio::abm::LocationType::Home, person), &home);
}

TEST(TestWorld, evolveStateTransition)
{
    using testing::Return;

    auto world     = mio::abm::World(6);
    auto location1 = world.add_location(mio::abm::LocationType::School);
    auto& p1       = world.add_person(location1, mio::abm::InfectionState::Carrier);
    auto& p2       = world.add_person(location1, mio::abm::InfectionState::Susceptible);
    auto location2 = world.add_location(mio::abm::LocationType::Work);
    auto& p3       = world.add_person(location2, mio::abm::InfectionState::Infected);
    p1.set_assigned_location(location1);
    p2.set_assigned_location(location1);
    p3.set_assigned_location(location2);

    //setup mock so only p2 transitions
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke)
        .Times(testing::AtLeast(3))
        .WillOnce(Return(0.51))
        .WillOnce(Return(0.04))
        .WillOnce(Return(0.6))
        .WillRepeatedly(Return(1.0));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    world.evolve(mio::abm::TimePoint(0), mio::abm::hours(1));

    EXPECT_EQ(p1.get_infection_state(), mio::abm::InfectionState::Carrier);
    EXPECT_EQ(p2.get_infection_state(), mio::abm::InfectionState::Exposed);
    EXPECT_EQ(p3.get_infection_state(), mio::abm::InfectionState::Infected);
}

TEST(TestWorld, evolveMigration)
{
    using testing::Return;

    {
        auto world     = mio::abm::World(6);
        auto home_id   = world.add_location(mio::abm::LocationType::Home);
        auto school_id = world.add_location(mio::abm::LocationType::School);
        auto work_id   = world.add_location(mio::abm::LocationType::Work);

        ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>>
            mock_uniform_dist;
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
            .Times(testing::Exactly(8))
            .WillOnce(testing::Return(0.8)) // draw random work group
            .WillOnce(testing::Return(0.8)) // draw random school group
            .WillOnce(testing::Return(0.8)) // draw random work hour
            .WillOnce(testing::Return(0.8)) // draw random school hour
            .WillOnce(testing::Return(0.8)) // draw random work group
            .WillOnce(testing::Return(0.8)) // draw random school group
            .WillOnce(testing::Return(0.8)) // draw random work hour
            .WillOnce(testing::Return(0.8)); // draw random school hour

        auto& p1 = world.add_person(home_id, mio::abm::InfectionState::Carrier, mio::AgeGroup(2));
        auto& p2 = world.add_person(home_id, mio::abm::InfectionState::Susceptible, mio::AgeGroup(1));

        p1.set_assigned_location(school_id);
        p2.set_assigned_location(school_id);
        p1.set_assigned_location(work_id);
        p2.set_assigned_location(work_id);
        p1.set_assigned_location(home_id);
        p2.set_assigned_location(home_id);

        auto& school = world.get_individualized_location(school_id);
        auto& work   = world.get_individualized_location(work_id);

        ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
            mock_exponential_dist;
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.)); //no state transitions

        world.evolve(mio::abm::TimePoint(0) + mio::abm::hours(8), mio::abm::hours(1));

        EXPECT_EQ(p1.get_location_id().type, mio::abm::LocationType::Work);
        EXPECT_EQ(p2.get_location_id().type, mio::abm::LocationType::School);
        EXPECT_EQ(school.get_population().get_last_value().sum(), 1);
        EXPECT_EQ(work.get_population().get_last_value().sum(), 1);
    }

    {
        auto world = mio::abm::World(6);
        world.use_migration_rules(false);

        auto home_id     = world.add_location(mio::abm::LocationType::Home);
        auto event_id    = world.add_location(mio::abm::LocationType::SocialEvent);
        auto work_id     = world.add_location(mio::abm::LocationType::Work);
        auto hospital_id = world.add_location(mio::abm::LocationType::Hospital);

        auto& p1 = world.add_person(home_id, mio::abm::InfectionState::Carrier, mio::AgeGroup(2));
        auto& p2 = world.add_person(home_id, mio::abm::InfectionState::Susceptible, mio::AgeGroup(1));
        auto& p3 = world.add_person(home_id, mio::abm::InfectionState::Infected_Severe, mio::AgeGroup(1));
        auto& p4 = world.add_person(hospital_id, mio::abm::InfectionState::Recovered_Infected, mio::AgeGroup(1));
        auto& p5 = world.add_person(home_id, mio::abm::InfectionState::Susceptible, mio::AgeGroup(2));
        p1.set_assigned_location(event_id);
        p2.set_assigned_location(event_id);
        p1.set_assigned_location(work_id);
        p2.set_assigned_location(work_id);
        p1.set_assigned_location(home_id);
        p2.set_assigned_location(home_id);
        p3.set_assigned_location(home_id);
        p4.set_assigned_location(home_id);
        p3.set_assigned_location(hospital_id);
        p5.set_assigned_location(event_id);
        p5.set_assigned_location(work_id);
        p5.set_assigned_location(home_id);

        mio::abm::TripList& data = world.get_trip_list();
        mio::abm::Trip trip1(p1.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), work_id, home_id);
        mio::abm::Trip trip2(p2.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), event_id, home_id);
        mio::abm::Trip trip3(p5.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), event_id, work_id);
        data.add_trip(trip1);
        data.add_trip(trip2);
        data.add_trip(trip3);

        ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
            mock_exponential_dist;
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.)); //no state transitions

        world.evolve(mio::abm::TimePoint(0) + mio::abm::hours(8), mio::abm::hours(2));

        auto& event    = world.get_individualized_location(event_id);
        auto& work     = world.get_individualized_location(work_id);
        auto& home     = world.get_individualized_location(home_id);
        auto& hospital = world.get_individualized_location(hospital_id);

        EXPECT_EQ(p1.get_location_id().type, mio::abm::LocationType::Work);
        EXPECT_EQ(p2.get_location_id().type, mio::abm::LocationType::SocialEvent);
        EXPECT_EQ(p3.get_location_id().type, mio::abm::LocationType::Hospital);
        EXPECT_EQ(p4.get_location_id().type, mio::abm::LocationType::Home);
        EXPECT_EQ(p5.get_location_id().type, mio::abm::LocationType::Home);
        EXPECT_EQ(event.get_population().get_last_value().sum(), 1);
        EXPECT_EQ(work.get_population().get_last_value().sum(), 1);
        EXPECT_EQ(home.get_population().get_last_value().sum(), 2);
        EXPECT_EQ(hospital.get_population().get_last_value().sum(), 1);
    }
}

TEST(TestWorldTestingCriteria, testAddingAndUpdatingAndRunningTestingSchemes)
{

    auto world   = mio::abm::World(6);
    auto home_id = world.add_location(mio::abm::LocationType::Home);
    auto work_id = world.add_location(mio::abm::LocationType::Work);
    auto person  = mio::abm::Person(home_id, mio::abm::InfectionState::Infected, mio::AgeGroup(2),
                                    world.get_global_infection_parameters());
    auto& home   = world.get_individualized_location(home_id);
    auto& work   = world.get_individualized_location(work_id);
    person.set_assigned_location(home);
    person.set_assigned_location(work);

    auto testing_criteria = mio::abm::TestingCriteria({}, {}, {});
    testing_criteria.add_infection_state(mio::abm::InfectionState::Infected);
    testing_criteria.add_infection_state(mio::abm::InfectionState::Carrier);
    testing_criteria.add_location_type(mio::abm::LocationType::Home);
    testing_criteria.add_location_type(mio::abm::LocationType::Work);

    const auto testing_frequency = mio::abm::days(1);
    const auto start_date        = mio::abm::TimePoint(20);
    const auto end_date          = mio::abm::TimePoint(60 * 60 * 24 * 3);
    const auto probability       = 1.0;
    const auto test_type         = mio::abm::PCRTest();

    auto testing_scheme =
        mio::abm::TestingScheme({testing_criteria}, testing_frequency, start_date, end_date, test_type, probability);

    world.get_testing_strategy().add_testing_scheme(testing_scheme);
    auto current_time = mio::abm::TimePoint(0);
    ASSERT_EQ(world.get_testing_strategy().run_strategy(person, work),
              true); // no active testing scheme -> person can enter
    current_time = mio::abm::TimePoint(30);
    world.get_testing_strategy().update_activity_status(current_time);
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.4));
    ASSERT_EQ(world.get_testing_strategy().run_strategy(person, work), false);

    world.get_testing_strategy().add_testing_scheme(testing_scheme); //doesn't get added because of == operator
    world.get_testing_strategy().remove_testing_scheme(testing_scheme);
    ASSERT_EQ(world.get_testing_strategy().run_strategy(person, work), true); // no more testing_schemes
}
