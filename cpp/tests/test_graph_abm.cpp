/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Julia Bicker, Martin J. Kuehn
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
#include "models/graph_abm/graph_world.h"
#include "models/graph_abm/graph_world.cpp"
#include "models/abm/world.h"
#include "models/abm/person.h"
#include "abm_helpers.h"
#include <gtest/gtest.h>

TEST(TestGraphWorld, test_find_location)
{
    mio::abm::GlobalInfectionParameters infection_params;
    auto world1 = mio::graph_abm::GraphWorld(infection_params, 0);
    auto world2 = mio::graph_abm::GraphWorld(infection_params, 1);
    auto home   = world1.add_location(mio::abm::LocationType::Home);
    auto work   = world2.add_location(mio::abm::LocationType::Work);
    auto person = mio::abm::Person(world1.get_individualized_location({0, mio::abm::LocationType::Home, 0}),
                                   mio::abm::AgeGroup::Age35to59, 0, 0);
    person.set_assigned_location(home);
    person.set_assigned_location(work);
    auto& found_home = world1.find_location(mio::abm::LocationType::Home, person);
    auto& found_work = world1.find_location(mio::abm::LocationType::Work, person);
    ASSERT_EQ(found_home, home);
    ASSERT_EQ(found_work, work);
}

TEST(TestGraphWorld, test_evolve_state_transition)
{
    using testing::Return;

    auto t  = mio::abm::TimePoint(0);
    auto dt = mio::abm::hours(1);

    auto params = mio::abm::GlobalInfectionParameters{};
    //setup so p1 and p3 don't transition
    params.get<mio::abm::IncubationPeriod>()[{mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
                                              mio::abm::VaccinationState::Unvaccinated}] = 2 * dt.days();
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{
        mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34, mio::abm::VaccinationState::Unvaccinated}] =
        2 * dt.days();
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{
        mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34, mio::abm::VaccinationState::Unvaccinated}] =
        2 * dt.days();
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
                                                      mio::abm::VaccinationState::Unvaccinated}] = 2 * dt.days();
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{
        mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34, mio::abm::VaccinationState::Unvaccinated}] =
        2 * dt.days();

    auto world     = mio::graph_abm::GraphWorld(params, 0);
    auto location1 = world.add_location(mio::abm::LocationType::School);
    auto& p1 =
        add_test_person(world, location1, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms);
    auto& p2 = add_test_person(world, location1, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::Susceptible);
    p1.set_assigned_location(location1);
    p2.set_assigned_location(location1);

    //setup mock so p2 becomes infected
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.0));

    world.evolve(t, dt);

    EXPECT_EQ(p1.get_infection_state(t + dt), mio::abm::InfectionState::InfectedSymptoms);
    EXPECT_EQ(p2.get_infection_state(t + dt), mio::abm::InfectionState::Exposed);
}

TEST(TestGraphWorld, test_evolve_migration)
{
    using testing::Return;

    { //test migration rules
        auto t  = mio::abm::TimePoint(0) + mio::abm::hours(6);
        auto dt = mio::abm::hours(1);

        auto params = mio::abm::GlobalInfectionParameters{};

        auto world1 = mio::graph_abm::GraphWorld(params, 0);
        auto world2 = mio::graph_abm::GraphWorld(params, 1);
        //home and school are in world1
        auto home_id   = world1.add_location(mio::abm::LocationType::Home);
        auto school_id = world1.add_location(mio::abm::LocationType::School);
        //work is in world2
        auto work_id = world2.add_location(mio::abm::LocationType::Work);

        ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>>
            mock_uniform_dist;
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
            .Times(testing::AtLeast(8))
            .WillOnce(testing::Return(0.0)) // draw random work group
            .WillOnce(testing::Return(0.0)) // draw random school group
            .WillOnce(testing::Return(0.0)) // draw random work hour
            .WillOnce(testing::Return(0.0)) // draw random school hour
            .WillOnce(testing::Return(0.0)) // draw random work group
            .WillOnce(testing::Return(0.0)) // draw random school group
            .WillOnce(testing::Return(0.0)) // draw random work hour
            .WillOnce(testing::Return(0.0)) // draw random school hour
            .WillRepeatedly(testing::Return(1.0));

        auto& p1 =
            add_test_person(world1, home_id, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::Susceptible, t);
        auto& p2 =
            add_test_person(world1, home_id, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::Susceptible, t);

        p1.set_assigned_location(home_id);
        p2.set_assigned_location(home_id);
        p1.set_assigned_location(work_id);
        p2.set_assigned_location(school_id);

        world1.evolve(t, dt);
        world2.evolve(t, dt);

        auto& school = world1.get_individualized_location(school_id);
        auto& work   = world1.get_individualized_location(work_id);

        EXPECT_EQ(p1.get_location(), work);
        EXPECT_EQ(p2.get_location(), school);
        EXPECT_EQ(school.get_number_persons(), 1);
        EXPECT_EQ(work.get_number_persons(), 1);

        EXPECT_EQ(world1.get_persons_to_migrate().size(), 1);
    }
    { //test trips
        auto t  = mio::abm::TimePoint(0) + mio::abm::hours(8);
        auto dt = mio::abm::hours(2);

        auto params = mio::abm::GlobalInfectionParameters{};

        //setup so p1-p5 don't transition
        params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{
            mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
            mio::abm::VaccinationState::Unvaccinated}]                                                  = 2 * dt.days();
        params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype,
                                                             mio::abm::AgeGroup::Age15to34,
                                                             mio::abm::VaccinationState::Unvaccinated}] = 2 * dt.days();
        params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
                                                  mio::abm::VaccinationState::Unvaccinated}]            = 2 * dt.days();
        params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
                                                   mio::abm::VaccinationState::Unvaccinated}]           = 2 * dt.days();

        auto world1 = mio::graph_abm::GraphWorld(params, 0);
        auto world2 = mio::graph_abm::GraphWorld(params, 1);
        world1.use_migration_rules(false);
        world2.use_migration_rules(false);
        //home is in world1
        auto home_id = world1.add_location(mio::abm::LocationType::Home);
        //other locations are in world2
        auto event_id    = world2.add_location(mio::abm::LocationType::SocialEvent);
        auto work_id     = world2.add_location(mio::abm::LocationType::Work);
        auto hospital_id = world2.add_location(mio::abm::LocationType::Hospital);

        //create persons
        auto& p1 =
            add_test_person(world1, home_id, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::Susceptible, t);
        auto& p2 =
            add_test_person(world1, home_id, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::Susceptible, t);
        auto& p3 =
            add_test_person(world1, home_id, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::InfectedSevere, t);
        auto& p4 =
            add_test_person(world1, home_id, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::Susceptible, t);
        auto& p5 =
            add_test_person(world1, home_id, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::Susceptible, t);

        p1.set_assigned_location(home_id);
        p1.set_assigned_location(work_id);
        p1.set_assigned_location(event_id);
        p2.set_assigned_location(home_id);
        p2.set_assigned_location(event_id);
        p3.set_assigned_location(home_id);
        p3.set_assigned_location(hospital_id);
        p4.set_assigned_location(home_id);
        p4.set_assigned_location(work_id);
        p4.set_assigned_location(event_id);
        p5.set_assigned_location(home_id);
        p5.set_assigned_location(work_id);
        p5.set_assigned_location(event_id);

        mio::abm::TripList& trip_list = world1.get_trip_list();
        mio::abm::Trip trip1(p1.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), work_id, home_id, {},
                             p1.get_world_id());
        mio::abm::Trip trip2(p2.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), event_id, home_id, {},
                             p2.get_world_id());
        mio::abm::Trip trip3(p4.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), event_id, home_id, {},
                             p4.get_world_id());
        mio::abm::Trip trip4(p5.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), event_id, work_id, {},
                             p5.get_world_id());
        trip_list.add_trip(trip1);
        trip_list.add_trip(trip2);
        trip_list.add_trip(trip3);
        trip_list.add_trip(trip4);

        mio::unused(trip1);
        mio::unused(trip_list);

        ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
            mock_exponential_dist;
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.)); //no state transitions

        world1.evolve(t, dt);
        world2.evolve(t, dt);

        auto& event    = world1.get_individualized_location(event_id);
        auto& work     = world1.get_individualized_location(work_id);
        auto& home     = world1.get_individualized_location(home_id);
        auto& hospital = world1.get_individualized_location(hospital_id);

        EXPECT_EQ(p1.get_location(), work);
        EXPECT_EQ(p2.get_location(), event);
        EXPECT_EQ(p3.get_location(), hospital);
        EXPECT_EQ(p4.get_location(), event);
        EXPECT_EQ(p5.get_location(), home);
        EXPECT_EQ(event.get_number_persons(), 2);
        EXPECT_EQ(work.get_number_persons(), 1);
        EXPECT_EQ(home.get_number_persons(), 1);
        EXPECT_EQ(hospital.get_number_persons(), 1);

        EXPECT_EQ(world1.get_persons_to_migrate().size(), 4);
    }
}

TEST(TestGraphWorld, test_add_existing_person)
{
    auto params = mio::abm::GlobalInfectionParameters{};
    auto world  = mio::graph_abm::GraphWorld(params, 0);

    auto loc        = mio::abm::Location(mio::abm::LocationType::Home, 0, 0);
    auto person_ptr = std::make_unique<mio::abm::Person>(loc, mio::abm::AgeGroup::Age35to59, 0, 0);
    world.add_existing_person(std::move(person_ptr));

    EXPECT_EQ(world.get_persons().size(), 1);
    auto& person = *world.get_person(0, 0);
    mio::unused(person);
    EXPECT_EQ(person->get_person_id(), 0U);
    EXPECT_EQ(person->get_world_id(), 0U);
    EXPECT_EQ(person->get_age(), mio::abm::AgeGroup::Age35to59);
}
