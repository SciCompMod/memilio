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
#include "abm_helpers.h"

TEST(TestWorld, init)
{
    auto world = mio::abm::World();

    EXPECT_EQ(world.get_locations().size(), 1);
    EXPECT_EQ(world.get_locations()[0].get_type(), mio::abm::LocationType::Cemetery);
    ASSERT_THAT(world.get_persons(), testing::ElementsAre());
}

TEST(TestWorld, addLocation)
{
    auto world      = mio::abm::World();
    auto school_id1 = world.add_location(mio::abm::LocationType::School);
    auto school_id2 = world.add_location(mio::abm::LocationType::School);
    auto work_id    = world.add_location(mio::abm::LocationType::Work);
    auto home_id    = world.add_location(mio::abm::LocationType::Home);

    ASSERT_EQ((int)school_id1.index, 1);
    ASSERT_EQ((int)school_id2.index, 2);

    auto& school1 = world.get_individualized_location(school_id1);
    auto& school2 = world.get_individualized_location(school_id2);
    auto& work    = world.get_individualized_location(work_id);
    auto& home    = world.get_individualized_location(home_id);

    size_t count_schools = 0;
    for (auto& loc : world.get_locations()) {
        if (loc.get_type() == mio::abm::LocationType::School) {
            count_schools++;
        }
    }
    ASSERT_EQ(count_schools, 2);

    ASSERT_EQ(world.get_locations()[1], school1);
    ASSERT_EQ(world.get_locations()[2], school2);
    ASSERT_EQ(world.get_locations()[3], work);
    ASSERT_EQ(world.get_locations()[4], home);
}

TEST(TestWorld, addPerson)
{
    auto world    = mio::abm::World();
    auto location = world.add_location(mio::abm::LocationType::School);

    auto& p1 = world.add_person(location, mio::abm::AgeGroup::Age15to34);
    auto& p2 = world.add_person(location, mio::abm::AgeGroup::Age35to59);

    ASSERT_EQ(world.get_persons().size(), 2);
    ASSERT_EQ(&world.get_persons()[0], &p1);
    ASSERT_EQ(&world.get_persons()[1], &p2);
}

TEST(TestWorld, getSubpopulationCombined)
{
    auto t       = mio::abm::TimePoint(0);
    auto world   = mio::abm::World();
    auto school1 = world.add_location(mio::abm::LocationType::School);
    auto school2 = world.add_location(mio::abm::LocationType::School);
    auto school3 = world.add_location(mio::abm::LocationType::School);
    add_test_person(world, school1, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedNoSymptoms);
    add_test_person(world, school1, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::Susceptible);
    add_test_person(world, school2, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::Susceptible);
    add_test_person(world, school2, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::Susceptible);
    add_test_person(world, school3, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedNoSymptoms);

    ASSERT_EQ(
        world.get_subpopulation_combined(t, mio::abm::InfectionState::Susceptible, mio::abm::LocationType::School), 3);
    ASSERT_EQ(world.get_subpopulation_combined(t, mio::abm::InfectionState::InfectedNoSymptoms,
                                               mio::abm::LocationType::School),
              2);
}

TEST(TestWorld, findLocation)
{
    auto world     = mio::abm::World();
    auto home_id   = world.add_location(mio::abm::LocationType::Home);
    auto school_id = world.add_location(mio::abm::LocationType::School);
    auto work_id   = world.add_location(mio::abm::LocationType::Work);
    auto& home     = world.get_individualized_location(home_id);
    auto& school   = world.get_individualized_location(school_id);
    auto& work     = world.get_individualized_location(work_id);
    auto person    = make_test_person(home);
    person.set_assigned_location(home);
    person.set_assigned_location(work);
    person.set_assigned_location(school);

    ASSERT_EQ(world.find_location(mio::abm::LocationType::Work, person), work);
    ASSERT_EQ(world.find_location(mio::abm::LocationType::School, person), school);
    ASSERT_EQ(world.find_location(mio::abm::LocationType::Home, person), home);
}

TEST(TestWorld, evolveStateTransition)
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
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{
        mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34, mio::abm::VaccinationState::Unvaccinated}] =
        2 * dt.days();
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
                                                      mio::abm::VaccinationState::Unvaccinated}] = 2 * dt.days();
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{
        mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34, mio::abm::VaccinationState::Unvaccinated}] =
        2 * dt.days();

    auto world     = mio::abm::World(params);
    auto location1 = world.add_location(mio::abm::LocationType::School);
    auto& p1 =
        add_test_person(world, location1, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedNoSymptoms);
    auto& p2 = add_test_person(world, location1, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::Susceptible);
    auto location2 = world.add_location(mio::abm::LocationType::Work);
    auto& p3 =
        add_test_person(world, location2, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms);
    p1.set_assigned_location(location1);
    p2.set_assigned_location(location1);
    p3.set_assigned_location(location2);

    //setup mock so p2 becomes infected
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.0));

    world.evolve(t, dt);

    EXPECT_EQ(p1.get_infection_state(t + dt), mio::abm::InfectionState::InfectedNoSymptoms);
    EXPECT_EQ(p2.get_infection_state(t + dt), mio::abm::InfectionState::Exposed);
    EXPECT_EQ(p3.get_infection_state(t + dt), mio::abm::InfectionState::InfectedSymptoms);
}

TEST(TestWorld, evolveMigration)
{
    using testing::Return;

    {
        auto t      = mio::abm::TimePoint(0) + mio::abm::hours(8);
        auto dt     = mio::abm::hours(1);
        auto params = mio::abm::GlobalInfectionParameters{};
        //setup so p1 doesn't transition
        params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{
            mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
            mio::abm::VaccinationState::Unvaccinated}] = 2 * dt.days();
        params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{
            mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
            mio::abm::VaccinationState::Unvaccinated}] = 2 * dt.days();

        auto world     = mio::abm::World(params);
        auto home_id   = world.add_location(mio::abm::LocationType::Home);
        auto school_id = world.add_location(mio::abm::LocationType::School);
        auto work_id   = world.add_location(mio::abm::LocationType::Work);

        ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>>
            mock_uniform_dist;
        EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
            .Times(testing::AtLeast(8))
            .WillOnce(testing::Return(0.8)) // draw random work group
            .WillOnce(testing::Return(0.8)) // draw random school group
            .WillOnce(testing::Return(0.8)) // draw random work hour
            .WillOnce(testing::Return(0.8)) // draw random school hour
            .WillOnce(testing::Return(0.8)) // draw random work group
            .WillOnce(testing::Return(0.8)) // draw random school group
            .WillOnce(testing::Return(0.8)) // draw random work hour
            .WillOnce(testing::Return(0.8)) // draw random school hour
            .WillRepeatedly(testing::Return(1.0));

        auto& p2 =
            add_test_person(world, home_id, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::Susceptible, t);
        auto& p1 = add_test_person(world, home_id, mio::abm::AgeGroup::Age15to34,
                                   mio::abm::InfectionState::InfectedNoSymptoms, t);

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

        world.evolve(t, dt);

        EXPECT_EQ(p1.get_location(), work);
        EXPECT_EQ(p2.get_location(), school);
        EXPECT_EQ(school.get_number_persons(), 1);
        EXPECT_EQ(work.get_number_persons(), 1);
    }

    {
        auto t      = mio::abm::TimePoint(0) + mio::abm::hours(8);
        auto dt     = mio::abm::hours(2);
        auto params = mio::abm::GlobalInfectionParameters{};
        //setup so p1-p5 don't transition
        params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{
            mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
            mio::abm::VaccinationState::Unvaccinated}]                                        = 2 * dt.days();
        params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{
            mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
            mio::abm::VaccinationState::Unvaccinated}]                                        = 2 * dt.days();
        params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
                                                  mio::abm::VaccinationState::Unvaccinated}]  = 2 * dt.days();
        params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age15to34,
                                                   mio::abm::VaccinationState::Unvaccinated}] = 2 * dt.days();

        auto world = mio::abm::World(params);
        world.use_migration_rules(false);

        auto home_id     = world.add_location(mio::abm::LocationType::Home);
        auto event_id    = world.add_location(mio::abm::LocationType::SocialEvent);
        auto work_id     = world.add_location(mio::abm::LocationType::Work);
        auto hospital_id = world.add_location(mio::abm::LocationType::Hospital);

        auto& p1 = add_test_person(world, home_id, mio::abm::AgeGroup::Age15to34,
                                   mio::abm::InfectionState::InfectedNoSymptoms, t);
        auto& p2 =
            add_test_person(world, home_id, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::Susceptible, t);
        auto& p3 =
            add_test_person(world, home_id, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::InfectedSevere, t);
        auto& p4 =
            add_test_person(world, hospital_id, mio::abm::AgeGroup::Age5to14, mio::abm::InfectionState::Recovered, t);
        auto& p5 =
            add_test_person(world, home_id, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::Susceptible, t);
        p1.set_assigned_location(event_id);
        p2.set_assigned_location(event_id);
        p1.set_assigned_location(work_id);
        p2.set_assigned_location(work_id);
        p1.set_assigned_location(home_id);
        p2.set_assigned_location(home_id);
        p3.set_assigned_location(home_id);
        p4.set_assigned_location(home_id);
        p3.set_assigned_location(hospital_id);
        p4.set_assigned_location(hospital_id);
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

        world.evolve(t, dt);

        auto& event    = world.get_individualized_location(event_id);
        auto& work     = world.get_individualized_location(work_id);
        auto& home     = world.get_individualized_location(home_id);
        auto& hospital = world.get_individualized_location(hospital_id);

        EXPECT_EQ(p1.get_location(), work);
        EXPECT_EQ(p2.get_location(), event);
        EXPECT_EQ(p3.get_location(), hospital);
        EXPECT_EQ(p4.get_location(), home);
        EXPECT_EQ(p5.get_location(), home);
        EXPECT_EQ(event.get_number_persons(), 1);
        EXPECT_EQ(work.get_number_persons(), 1);
        EXPECT_EQ(home.get_number_persons(), 2);
        EXPECT_EQ(hospital.get_number_persons(), 1);
    }

    // Test that a dead person cannot make a movement
    {
        auto t      = mio::abm::TimePoint(0);
        auto dt     = mio::abm::days(1);
        auto params = mio::abm::GlobalInfectionParameters{};
        // Time to go from severe to critical infection is 1 day (dt).
        params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age60to79,
                                                  mio::abm::VaccinationState::Unvaccinated}] = 1;
        // Time to go from critical infection to dead state is 1/2 day (0.5 * dt).
        params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age60to79,
                                                mio::abm::VaccinationState::Unvaccinated}] = 0.5 * 1;

        auto world       = mio::abm::World(params);
        auto home_id     = world.add_location(mio::abm::LocationType::Home);
        auto work_id     = world.add_location(mio::abm::LocationType::Work);
        auto icu_id      = world.add_location(mio::abm::LocationType::ICU);
        auto hospital_id = world.add_location(mio::abm::LocationType::Hospital);
        // Create a person that is dead at time t
        auto& p_dead = add_test_person(world, icu_id, mio::abm::AgeGroup::Age60to79, mio::abm::InfectionState::Dead, t);
        // Create a person that is severe at hospital and will be dead at time t + dt
        auto& p_severe =
            add_test_person(world, hospital_id, mio::abm::AgeGroup::Age60to79, mio::abm::InfectionState::Dead, t + dt);
        p_dead.set_assigned_location(icu_id);
        p_dead.set_assigned_location(work_id);
        p_dead.set_assigned_location(home_id);
        p_severe.set_assigned_location(hospital_id);
        p_severe.set_assigned_location(icu_id);
        p_severe.set_assigned_location(home_id);

        // Add trip to see if a dead person can move outside of cemetery by scheduled
        mio::abm::TripList& trip_list = world.get_trip_list();
        mio::abm::Trip trip1(p_dead.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(2), work_id, home_id);
        mio::abm::Trip trip2(p_dead.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(3), home_id, icu_id);
        mio::abm::Trip trip3(p_severe.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(3), home_id, icu_id);
        trip_list.add_trip(trip1);
        trip_list.add_trip(trip2);
        trip_list.add_trip(trip3);

        // Check the dead person got burried and the severely infected person starts in Hospital
        world.evolve(t, dt);
        EXPECT_EQ(p_dead.get_location().get_type(), mio::abm::LocationType::Cemetery);
        EXPECT_EQ(p_severe.get_infection_state(t), mio::abm::InfectionState::InfectedSevere);
        EXPECT_EQ(p_severe.get_location().get_type(), mio::abm::LocationType::Hospital);

        // Check the dead person is still in Cemetery and the severely infected person dies and got burried
        world.evolve(t + dt, dt);
        EXPECT_EQ(p_dead.get_location().get_type(), mio::abm::LocationType::Cemetery);
        EXPECT_EQ(p_severe.get_infection_state(t + dt), mio::abm::InfectionState::Dead);
        EXPECT_EQ(p_severe.get_location().get_type(), mio::abm::LocationType::Cemetery);
    }
}

TEST(TestWorldTestingCriteria, testAddingAndUpdatingAndRunningTestingSchemes)
{
    mio::abm::GlobalInfectionParameters params;
    // make sure the infected person stay in Infected long enough
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant(0), mio::abm::AgeGroup::Age15to34,
                                                         mio::abm::VaccinationState::Unvaccinated}] = 100;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant(0), mio::abm::AgeGroup::Age15to34,
                                                      mio::abm::VaccinationState::Unvaccinated}]    = 100;

    auto world        = mio::abm::World(params);
    auto home_id      = world.add_location(mio::abm::LocationType::Home);
    auto work_id      = world.add_location(mio::abm::LocationType::Work);
    auto& home        = world.get_individualized_location(home_id);
    auto& work        = world.get_individualized_location(work_id);
    auto current_time = mio::abm::TimePoint(0);
    auto person       = add_test_person(world, home_id, mio::abm::AgeGroup::Age15to34,
                                        mio::abm::InfectionState::InfectedSymptoms, current_time);
    person.set_assigned_location(home);
    person.set_assigned_location(work);

    auto testing_criteria = mio::abm::TestingCriteria({}, {}, {});
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedSymptoms);
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedNoSymptoms);
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
    ASSERT_EQ(world.get_testing_strategy().run_strategy(person, work, current_time),
              true); // no active testing scheme -> person can enter
    current_time = mio::abm::TimePoint(30);
    world.get_testing_strategy().update_activity_status(current_time);
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.4));
    ASSERT_EQ(world.get_testing_strategy().run_strategy(person, work, current_time), false);

    world.get_testing_strategy().add_testing_scheme(testing_scheme); //doesn't get added because of == operator
    world.get_testing_strategy().remove_testing_scheme(testing_scheme);
    ASSERT_EQ(world.get_testing_strategy().run_strategy(person, work, current_time), true); // no more testing_schemes
}
