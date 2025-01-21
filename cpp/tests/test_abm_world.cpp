/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, David Kerkmann, Sascha Korf, Martin J. Kuehn, Khoa Nguyen, Carlotta Gerstein
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
#include "abm/person.h"
#include "abm_helpers.h"
#include "memilio/utils/random_number_generator.h"

TEST(TestWorld, init)
{
    auto world = mio::abm::World(num_age_groups);

    EXPECT_EQ(world.get_locations().size(), 1);
    EXPECT_EQ(world.get_locations()[0].get_type(), mio::abm::LocationType::Cemetery);
    ASSERT_THAT(world.get_persons(), testing::ElementsAre());
}

TEST(TestWorld, addLocation)
{
    auto world      = mio::abm::World(num_age_groups);
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
    auto world    = mio::abm::World(num_age_groups);
    auto location = world.add_location(mio::abm::LocationType::School);

    auto& p1 = world.add_person(location, age_group_15_to_34);
    auto& p2 = world.add_person(location, age_group_35_to_59);

    ASSERT_EQ(world.get_persons().size(), 2);
    ASSERT_EQ(&world.get_persons()[0], &p1);
    ASSERT_EQ(&world.get_persons()[1], &p2);
}

TEST(TestWorld, getSubpopulationCombined)
{
    auto t       = mio::abm::TimePoint(0);
    auto world   = mio::abm::World(num_age_groups);
    auto school1 = world.add_location(mio::abm::LocationType::School);
    auto school2 = world.add_location(mio::abm::LocationType::School);
    auto school3 = world.add_location(mio::abm::LocationType::School);
    auto home1   = world.add_location(mio::abm::LocationType::Home);
    add_test_person(world, school1, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);
    add_test_person(world, school1, age_group_15_to_34, mio::abm::InfectionState::Susceptible);
    add_test_person(world, school2, age_group_15_to_34, mio::abm::InfectionState::Susceptible);
    add_test_person(world, school2, age_group_15_to_34, mio::abm::InfectionState::Susceptible);
    add_test_person(world, school3, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);
    add_test_person(world, home1, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);

    ASSERT_EQ(world.get_subpopulation_combined_per_location_type(t, mio::abm::InfectionState::Susceptible,
                                                                 mio::abm::LocationType::School),
              3);
    ASSERT_EQ(world.get_subpopulation_combined_per_location_type(t, mio::abm::InfectionState::InfectedNoSymptoms,
                                                                 mio::abm::LocationType::School),
              2);
    ASSERT_EQ(world.get_subpopulation_combined(t, mio::abm::InfectionState::InfectedNoSymptoms), 3);
}

TEST(TestWorld, findLocation)
{
    auto world     = mio::abm::World(num_age_groups);
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

    auto&& world_test = std::as_const(world);
    ASSERT_EQ(world_test.find_location(mio::abm::LocationType::Work, person), work);
    ASSERT_EQ(world_test.find_location(mio::abm::LocationType::School, person), school);
    ASSERT_EQ(world_test.find_location(mio::abm::LocationType::Home, person), home);
}

TEST(TestWorld, evolveStateTransition)
{
    using testing::Return;

    auto t     = mio::abm::TimePoint(0);
    auto dt    = mio::abm::hours(1);
    auto world = mio::abm::World(num_age_groups);

    //setup so p1 and p3 don't transition
    world.parameters.get<mio::abm::IncubationPeriod>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        2 * dt.days();
    world.parameters
        .get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        2 * dt.days();
    world.parameters
        .get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        2 * dt.days();
    world.parameters.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        2 * dt.days();
    world.parameters
        .get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        2 * dt.days();

    auto location1 = world.add_location(mio::abm::LocationType::School);
    auto& p1 = add_test_person(world, location1, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);
    auto& p2 = add_test_person(world, location1, age_group_15_to_34, mio::abm::InfectionState::Susceptible);
    auto location2 = world.add_location(mio::abm::LocationType::Work);
    auto& p3       = add_test_person(world, location2, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms);
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
        auto t     = mio::abm::TimePoint(0) + mio::abm::hours(8);
        auto dt    = mio::abm::hours(1);
        auto world = mio::abm::World(num_age_groups);
        //setup so p1 doesn't do transition
        world.parameters
            .get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
            2 * dt.days();
        world.parameters
            .get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
            2 * dt.days();

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

        auto& p2 = add_test_person(world, home_id, age_group_5_to_14, mio::abm::InfectionState::Susceptible, t);
        auto& p1 = add_test_person(world, home_id, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms, t);

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
        auto t     = mio::abm::TimePoint(0) + mio::abm::hours(8);
        auto dt    = mio::abm::hours(2);
        auto world = mio::abm::World(num_age_groups);
        //setup so p1-p5 don't do transition
        world.parameters
            .get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
            2 * dt.days();
        world.parameters
            .get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
            2 * dt.days();
        world.parameters.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
            2 * dt.days();
        world.parameters.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
            2 * dt.days();

        auto home_id     = world.add_location(mio::abm::LocationType::Home);
        auto event_id    = world.add_location(mio::abm::LocationType::SocialEvent);
        auto work_id     = world.add_location(mio::abm::LocationType::Work);
        auto hospital_id = world.add_location(mio::abm::LocationType::Hospital);

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
            .WillRepeatedly(testing::Return(1.0)); // this forces p1 and p3 to recover

        auto& p1 = add_test_person(world, home_id, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms, t);
        auto& p2 = add_test_person(world, home_id, age_group_5_to_14, mio::abm::InfectionState::Susceptible, t);
        auto& p3 = add_test_person(world, home_id, age_group_5_to_14, mio::abm::InfectionState::InfectedSevere, t);
        auto& p4 = add_test_person(world, hospital_id, age_group_5_to_14, mio::abm::InfectionState::Recovered, t);
        auto& p5 = add_test_person(world, home_id, age_group_15_to_34, mio::abm::InfectionState::Susceptible, t);
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
        mio::abm::Trip trip3(p5.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), event_id, home_id);
        data.add_trip(trip1);
        data.add_trip(trip2);
        data.add_trip(trip3);

        data.use_weekday_trips_on_weekend();

        ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
            mock_exponential_dist;
        EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.)); //no infections

        world.evolve(t, dt);

        auto& event    = world.get_individualized_location(event_id);
        auto& work     = world.get_individualized_location(work_id);
        auto& home     = world.get_individualized_location(home_id);
        auto& hospital = world.get_individualized_location(hospital_id);

        EXPECT_EQ(p1.get_location(), work);
        EXPECT_EQ(p2.get_location(), event);
        EXPECT_EQ(p3.get_location(), hospital);
        EXPECT_EQ(p4.get_location(), home);
        EXPECT_EQ(p5.get_location(), event);
        EXPECT_EQ(event.get_number_persons(), 2);
        EXPECT_EQ(work.get_number_persons(), 1);
        EXPECT_EQ(home.get_number_persons(), 1);
        EXPECT_EQ(hospital.get_number_persons(), 1);

        p1.migrate_to(home);
        p2.migrate_to(home);
        p5.migrate_to(home);

        t = mio::abm::TimePoint(0) + mio::abm::days(6) + mio::abm::hours(8);
        world.get_trip_list().reset_index();

        world.evolve(t, dt);

        EXPECT_EQ(p1.get_location(), work);
        EXPECT_EQ(p2.get_location(), event);
        EXPECT_EQ(p3.get_location(), home);
        EXPECT_EQ(p4.get_location(), home);
        EXPECT_EQ(p5.get_location(), event);
        EXPECT_EQ(event.get_number_persons(), 2);
        EXPECT_EQ(work.get_number_persons(), 1);
        EXPECT_EQ(home.get_number_persons(), 2);

        bool weekend = true;
        mio::abm::Trip tripweekend1(p1.get_person_id(),
                                    mio::abm::TimePoint(0) + mio::abm::days(6) + mio::abm::hours(10), event_id);
        mio::abm::Trip tripweekend2(p2.get_person_id(),
                                    mio::abm::TimePoint(0) + mio::abm::days(6) + mio::abm::hours(10), home_id);
        mio::abm::Trip tripweekend3(p5.get_person_id(),
                                    mio::abm::TimePoint(0) + mio::abm::days(6) + mio::abm::hours(10), work_id);
        data.add_trip(tripweekend1, weekend);
        data.add_trip(tripweekend2, weekend);
        data.add_trip(tripweekend3, weekend);

        t += mio::abm::hours(1);

        world.evolve(t, dt);

        EXPECT_EQ(p1.get_location(), event);
        EXPECT_EQ(p2.get_location(), home);
        EXPECT_EQ(p3.get_location(), home);
        EXPECT_EQ(p4.get_location(), home);
        EXPECT_EQ(p5.get_location(), work);
        EXPECT_EQ(event.get_number_persons(), 1);
        EXPECT_EQ(work.get_number_persons(), 1);
        EXPECT_EQ(home.get_number_persons(), 3);
    }

    // Test that a dead person cannot make a movement
    {
        auto t     = mio::abm::TimePoint(0);
        auto dt    = mio::abm::days(1);
        auto world = mio::abm::World(num_age_groups);

        // Time to go from severe to critical infection is 1 day (dt).
        world.parameters.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] =
            0.5;
        // Time to go from critical infection to dead state is 1/2 day (0.5 * dt).
        world.parameters.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.5;

        auto home_id     = world.add_location(mio::abm::LocationType::Home);
        auto work_id     = world.add_location(mio::abm::LocationType::Work);
        auto icu_id      = world.add_location(mio::abm::LocationType::ICU);
        auto hospital_id = world.add_location(mio::abm::LocationType::Hospital);
        // Create a person that is dead at time t
        auto& p_dead = add_test_person(world, icu_id, age_group_60_to_79, mio::abm::InfectionState::Dead, t);
        // Create a person that is severe at hospital and will be dead at time t + dt
        auto& p_severe =
            add_test_person(world, hospital_id, age_group_60_to_79, mio::abm::InfectionState::Dead, t + dt);
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
    auto rng = mio::RandomNumberGenerator();

    auto world = mio::abm::World(num_age_groups);
    // make sure the infected person stay in Infected long enough
    world.parameters.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant(0), age_group_15_to_34}] =
        100;
    world.parameters.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant(0), age_group_15_to_34}] = 100;

    auto home_id      = world.add_location(mio::abm::LocationType::Home);
    auto work_id      = world.add_location(mio::abm::LocationType::Work);
    auto& home        = world.get_individualized_location(home_id);
    auto& work        = world.get_individualized_location(work_id);
    auto current_time = mio::abm::TimePoint(0);
    auto person =
        add_test_person(world, home_id, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms, current_time);
    auto rng_person = mio::abm::Person::RandomNumberGenerator(rng, person);
    person.set_assigned_location(home);
    person.set_assigned_location(work);

    auto testing_criteria = mio::abm::TestingCriteria();
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedSymptoms);
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedNoSymptoms);

    const auto testing_frequency = mio::abm::days(1);
    const auto start_date        = mio::abm::TimePoint(20);
    const auto end_date          = mio::abm::TimePoint(60 * 60 * 24 * 3);
    const auto probability       = 1.0;
    const auto test_type         = mio::abm::PCRTest();

    auto testing_scheme =
        mio::abm::TestingScheme(testing_criteria, testing_frequency, start_date, end_date, test_type, probability);

    world.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme);
    ASSERT_EQ(world.get_testing_strategy().run_strategy(rng_person, person, work, current_time),
              true); // no active testing scheme -> person can enter
    current_time = mio::abm::TimePoint(30);
    world.get_testing_strategy().update_activity_status(current_time);
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.7))
        .WillOnce(testing::Return(0.4));
    ASSERT_EQ(world.get_testing_strategy().run_strategy(rng_person, person, work, current_time), false);

    world.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work,
                                                    testing_scheme); //doesn't get added because of == operator
    world.get_testing_strategy().remove_testing_scheme(mio::abm::LocationType::Work, testing_scheme);
    ASSERT_EQ(world.get_testing_strategy().run_strategy(rng_person, person, work, current_time),
              true); // no more testing_schemes
}

TEST(TestWorld, checkParameterConstraints)
{
    mio::set_log_level(mio::LogLevel::critical); //errors inevitable as these are wanted
    auto world  = mio::abm::World(num_age_groups);
    auto params = world.parameters;

    params.get<mio::abm::IncubationPeriod>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]              = 1.;
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]  = 2.;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = 3.;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 4.;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]      = 5.;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]              = 6.;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]             = 7.;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]                = 8.;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]           = 9.;
    params.get<mio::abm::RecoveredToSusceptible>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]        = 10.;
    params.get<mio::abm::DetectInfection>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]               = 0.3;
    params.get<mio::abm::GotoWorkTimeMinimum>()[age_group_35_to_59]       = mio::abm::hours(4);
    params.get<mio::abm::GotoWorkTimeMaximum>()[age_group_35_to_59]       = mio::abm::hours(8);
    params.get<mio::abm::GotoSchoolTimeMinimum>()[age_group_0_to_4]       = mio::abm::hours(3);
    params.get<mio::abm::GotoSchoolTimeMaximum>()[age_group_0_to_4]       = mio::abm::hours(6);
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::Community] = 0.5;
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::FFP2]      = 0.6;
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::Surgical]  = 0.7;
    params.get<mio::abm::LockdownDate>()                                  = mio::abm::TimePoint(0);
    ASSERT_EQ(params.check_constraints(), false);

    params.get<mio::abm::IncubationPeriod>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = -1.;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::IncubationPeriod>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]             = 1.;
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = -2.;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::InfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]  = 2.;
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = -3.;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::InfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = 3.;
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = -4.;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = 4.;
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]    = -5.;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::InfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = 5.;
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]         = -6.;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::SevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]  = 6.;
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = -7.;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::SevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = 7.;
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]    = -8.;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::CriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]      = 8.;
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = -9.;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::CriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]    = 9.;
    params.get<mio::abm::RecoveredToSusceptible>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = -10.;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::RecoveredToSusceptible>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = 10.;
    params.get<mio::abm::DetectInfection>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]        = 1.1;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::DetectInfection>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = 0.3;

    params.get<mio::abm::GotoWorkTimeMinimum>()[age_group_35_to_59] = mio::abm::hours(30);
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::GotoWorkTimeMinimum>()[age_group_35_to_59] = mio::abm::hours(4);
    params.get<mio::abm::GotoWorkTimeMaximum>()[age_group_35_to_59] = mio::abm::hours(30);
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::GotoWorkTimeMaximum>()[age_group_35_to_59] = mio::abm::hours(8);
    params.get<mio::abm::GotoSchoolTimeMinimum>()[age_group_0_to_4] = mio::abm::hours(30);
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::GotoSchoolTimeMinimum>()[age_group_0_to_4] = mio::abm::hours(3);
    params.get<mio::abm::GotoSchoolTimeMaximum>()[age_group_0_to_4] = mio::abm::hours(30);
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::GotoSchoolTimeMaximum>()[age_group_0_to_4] = mio::abm::hours(6);

    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::Community] = 1.2;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::Community] = 0.5;
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::FFP2]      = 1.2;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::FFP2]     = 0.6;
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::Surgical] = 1.2;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::Surgical] = 0.7;

    params.get<mio::abm::LockdownDate>() = mio::abm::TimePoint(-2);
    ASSERT_EQ(params.check_constraints(), true);
}

TEST(TestWorld, copyWorld)
{
    auto world = mio::abm::World(num_age_groups);
    auto rng   = mio::RandomNumberGenerator();

    world.parameters.get<mio::abm::IncubationPeriod>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = 4.;
    world.use_migration_rules(false);

    auto school_id1 = world.add_location(mio::abm::LocationType::School);
    auto school_id2 = world.add_location(mio::abm::LocationType::School);
    auto work_id    = world.add_location(mio::abm::LocationType::Work);
    auto home_id    = world.add_location(mio::abm::LocationType::Home);

    auto& school1 = world.get_individualized_location(school_id1);
    school1.set_required_mask(mio::abm::MaskType::Surgical);
    school1.set_npi_active(true);
    auto& school2 = world.get_individualized_location(school_id2);
    school2.set_required_mask(mio::abm::MaskType::FFP2);
    auto& work = world.get_individualized_location(work_id);
    auto& home = world.get_individualized_location(home_id);

    auto& p1    = world.add_person(school_id1, age_group_0_to_4);
    auto rng_p1 = mio::abm::Person::RandomNumberGenerator(rng, p1);
    p1.add_new_infection(mio::abm::Infection(rng_p1, mio::abm::VirusVariant::Wildtype, p1.get_age(), world.parameters,
                                             mio::abm::TimePoint(0)));
    auto& p2 = world.add_person(school_id2, age_group_15_to_34);
    p2.set_mask_preferences(std::vector<double>(15, 0.2));

    mio::abm::TripList& trip_data = world.get_trip_list();
    mio::abm::Trip trip1(p1.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(8), school_id1, home_id);
    mio::abm::Trip trip2(p2.get_person_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), work_id, home_id);
    trip_data.add_trip(trip1);
    trip_data.add_trip(trip2);

    auto infection_params =
        world.parameters.get<mio::abm::IncubationPeriod>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]
            .value();

    auto copied_world = mio::abm::World(world);
    auto copied_infection_params =
        copied_world.parameters.get<mio::abm::IncubationPeriod>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]
            .value();

    // Assert the parameters, trips, locations and persons of copied world are logically equal to that of original world
    ASSERT_EQ(copied_infection_params, infection_params);
    ASSERT_EQ(copied_world.use_migration_rules(), world.use_migration_rules());

    mio::abm::TripList& copied_trip_data = copied_world.get_trip_list();
    ASSERT_EQ(copied_trip_data.num_trips(), trip_data.num_trips());
    ASSERT_EQ(copied_trip_data.get_next_trip(false).person_id, trip_data.get_next_trip(false).person_id);
    ASSERT_EQ(copied_trip_data.get_next_trip(false).migration_destination,
              trip_data.get_next_trip(false).migration_destination);
    ASSERT_EQ(copied_trip_data.get_next_trip(false).migration_origin, trip_data.get_next_trip(false).migration_origin);
    copied_trip_data.increase_index();
    trip_data.increase_index();
    ASSERT_EQ(copied_trip_data.get_next_trip(false).person_id, trip_data.get_next_trip(false).person_id);
    ASSERT_EQ(copied_trip_data.get_next_trip(false).migration_destination,
              trip_data.get_next_trip(false).migration_destination);
    ASSERT_EQ(copied_trip_data.get_next_trip(false).migration_origin, trip_data.get_next_trip(false).migration_origin);

    ASSERT_EQ(copied_world.get_locations().size(), world.get_locations().size());
    ASSERT_EQ(copied_world.get_locations()[1].get_index(), world.get_locations()[1].get_index());
    ASSERT_EQ(copied_world.get_locations()[2].get_index(), world.get_locations()[2].get_index());
    ASSERT_EQ(copied_world.get_locations()[3].get_index(), world.get_locations()[3].get_index());
    ASSERT_EQ(copied_world.get_locations()[4].get_index(), world.get_locations()[4].get_index());
    ASSERT_EQ(copied_world.get_locations()[1].get_number_persons(), world.get_locations()[1].get_number_persons());
    ASSERT_EQ(copied_world.get_locations()[2].get_number_persons(), world.get_locations()[2].get_number_persons());
    ASSERT_EQ(copied_world.get_locations()[3].get_number_persons(), world.get_locations()[3].get_number_persons());
    ASSERT_EQ(copied_world.get_locations()[4].get_number_persons(), world.get_locations()[4].get_number_persons());
    ASSERT_EQ(copied_world.get_locations()[1].get_npi_active(), world.get_locations()[1].get_npi_active());
    ASSERT_EQ(copied_world.get_locations()[2].get_npi_active(), world.get_locations()[2].get_npi_active());
    ASSERT_EQ(copied_world.get_locations()[3].get_npi_active(), world.get_locations()[3].get_npi_active());
    ASSERT_EQ(copied_world.get_locations()[4].get_npi_active(), world.get_locations()[4].get_npi_active());
    ASSERT_EQ(copied_world.get_locations()[1].get_required_mask(), world.get_locations()[1].get_required_mask());
    ASSERT_EQ(copied_world.get_locations()[2].get_required_mask(), world.get_locations()[2].get_required_mask());
    ASSERT_EQ(
        copied_world.get_locations()[1].get_subpopulation(mio::abm::TimePoint(0), mio::abm::InfectionState::Exposed),
        world.get_locations()[1].get_subpopulation(mio::abm::TimePoint(0), mio::abm::InfectionState::Exposed));
    ASSERT_EQ(
        copied_world.get_locations()[1].get_subpopulation(mio::abm::TimePoint(0),
                                                          mio::abm::InfectionState::Susceptible),
        world.get_locations()[1].get_subpopulation(mio::abm::TimePoint(0), mio::abm::InfectionState::Susceptible));
    ASSERT_EQ(
        copied_world.get_locations()[2].get_subpopulation(mio::abm::TimePoint(0), mio::abm::InfectionState::Exposed),
        world.get_locations()[2].get_subpopulation(mio::abm::TimePoint(0), mio::abm::InfectionState::Exposed));
    ASSERT_EQ(
        copied_world.get_locations()[2].get_subpopulation(mio::abm::TimePoint(0),
                                                          mio::abm::InfectionState::Susceptible),
        world.get_locations()[2].get_subpopulation(mio::abm::TimePoint(0), mio::abm::InfectionState::Susceptible));
    ASSERT_EQ(
        copied_world.get_locations()[3].get_subpopulation(mio::abm::TimePoint(0), mio::abm::InfectionState::Exposed),
        world.get_locations()[3].get_subpopulation(mio::abm::TimePoint(0), mio::abm::InfectionState::Exposed));
    ASSERT_EQ(
        copied_world.get_locations()[4].get_subpopulation(mio::abm::TimePoint(0), mio::abm::InfectionState::Exposed),
        world.get_locations()[4].get_subpopulation(mio::abm::TimePoint(0), mio::abm::InfectionState::Exposed));
    ASSERT_EQ(copied_world.get_locations()[1].get_cells().size(), world.get_locations()[1].get_cells().size());
    ASSERT_EQ(copied_world.get_locations()[2].get_cells().size(), world.get_locations()[2].get_cells().size());
    ASSERT_EQ(copied_world.get_locations()[3].get_cells().size(), world.get_locations()[2].get_cells().size());
    ASSERT_EQ(copied_world.get_locations()[4].get_cells().size(), world.get_locations()[2].get_cells().size());
    ASSERT_EQ(copied_world.get_locations()[1].get_cells()[0].m_persons.size(),
              world.get_locations()[1].get_cells()[0].m_persons.size());
    ASSERT_EQ(copied_world.get_locations()[2].get_cells()[0].m_persons.size(),
              world.get_locations()[2].get_cells()[0].m_persons.size());
    ASSERT_EQ(copied_world.get_locations()[3].get_cells()[0].m_persons.size(),
              world.get_locations()[3].get_cells()[0].m_persons.size());
    ASSERT_EQ(copied_world.get_locations()[4].get_cells()[0].m_persons.size(),
              world.get_locations()[4].get_cells()[0].m_persons.size());
    ASSERT_EQ(copied_world.get_locations()[1].get_cells()[0].m_persons[0],
              world.get_locations()[1].get_cells()[0].m_persons[0]);
    ASSERT_EQ(copied_world.get_locations()[2].get_cells()[0].m_persons[0],
              world.get_locations()[2].get_cells()[0].m_persons[0]);

    ASSERT_EQ(copied_world.get_persons().size(), world.get_persons().size());
    ASSERT_EQ(copied_world.get_persons()[0].get_location().get_index(),
              world.get_persons()[0].get_location().get_index());
    ASSERT_EQ(copied_world.get_persons()[1].get_location().get_index(),
              world.get_persons()[1].get_location().get_index());
    ASSERT_EQ(copied_world.get_persons()[0].get_location().get_type(),
              world.get_persons()[0].get_location().get_type());
    ASSERT_EQ(copied_world.get_persons()[1].get_location().get_type(),
              world.get_persons()[1].get_location().get_type());
    ASSERT_EQ(copied_world.get_persons()[0].get_infection().get_infection_state(mio::abm::TimePoint(0)),
              world.get_persons()[0].get_infection().get_infection_state(mio::abm::TimePoint(0)));
    ASSERT_EQ(copied_world.get_persons()[0].get_mask_compliance(mio::abm::LocationType::Home),
              world.get_persons()[0].get_mask_compliance(mio::abm::LocationType::Home));
    ASSERT_EQ(copied_world.get_persons()[0].get_mask_compliance(mio::abm::LocationType::Work),
              world.get_persons()[0].get_mask_compliance(mio::abm::LocationType::Work));
    ASSERT_EQ(copied_world.get_persons()[1].get_mask_compliance(mio::abm::LocationType::Home),
              world.get_persons()[1].get_mask_compliance(mio::abm::LocationType::Home));
    ASSERT_EQ(copied_world.get_persons()[1].get_mask_compliance(mio::abm::LocationType::Work),
              world.get_persons()[1].get_mask_compliance(mio::abm::LocationType::Work));

    // Assert the parameters, trips, locations, persons and their member variables of copied world are stored in different address of original world
    ASSERT_NE(&(copied_world.parameters), &world.parameters);
    ASSERT_NE(&(copied_world.get_trip_list()), &trip_data);

    ASSERT_NE(&copied_world.get_locations()[1], &world.get_locations()[1]);
    ASSERT_NE(&copied_world.get_locations()[2], &world.get_locations()[2]);
    ASSERT_NE(&copied_world.get_locations()[3], &world.get_locations()[3]);
    ASSERT_NE(&copied_world.get_locations()[4], &world.get_locations()[4]);
    ASSERT_NE(&copied_world.get_locations()[1].get_cells(), &world.get_locations()[1].get_cells());
    ASSERT_NE(&copied_world.get_locations()[2].get_cells(), &world.get_locations()[2].get_cells());
    ASSERT_NE(&copied_world.get_locations()[3].get_cells(), &world.get_locations()[3].get_cells());
    ASSERT_NE(&copied_world.get_locations()[4].get_cells(), &world.get_locations()[4].get_cells());
    ASSERT_NE(&(copied_world.get_locations()[1].get_cells()[0]), &(world.get_locations()[1].get_cells()[0]));
    ASSERT_NE(&(copied_world.get_locations()[2].get_cells()[0]), &(world.get_locations()[2].get_cells()[0]));
    ASSERT_NE(&(copied_world.get_locations()[1].get_cells()[0].m_persons[0]),
              &(world.get_locations()[1].get_cells()[0].m_persons[0]));
    ASSERT_NE(&(copied_world.get_locations()[2].get_cells()[0].m_persons[0]),
              &(world.get_locations()[2].get_cells()[0].m_persons[0]));

    ASSERT_NE(&copied_world.get_persons()[0], &world.get_persons()[0]);
    ASSERT_NE(&copied_world.get_persons()[1], &world.get_persons()[1]);
    ASSERT_NE(&(copied_world.get_persons()[0].get_location()), &world.get_persons()[0].get_location());
    ASSERT_NE(&(copied_world.get_persons()[1].get_location()), &world.get_persons()[1].get_location());
    ASSERT_NE(&(copied_world.get_locations()[1]), &(world.get_locations()[1]));
    ASSERT_NE(&(copied_world.get_locations()[2]), &(world.get_locations()[2]));
    ASSERT_NE(&(copied_world.get_persons()[0].get_assigned_locations()),
              &world.get_persons()[0].get_assigned_locations());
    ASSERT_NE(&(copied_world.get_persons()[1].get_assigned_locations()),
              &world.get_persons()[1].get_assigned_locations());
    ASSERT_NE(&(copied_world.get_persons()[0].get_infection()), &world.get_persons()[0].get_infection());
    ASSERT_NE(&(copied_world.get_persons()[0].get_mask()), &world.get_persons()[0].get_mask());
    ASSERT_NE(&(copied_world.get_persons()[1].get_mask()), &world.get_persons()[1].get_mask());
    ASSERT_NE(&(copied_world.get_persons()[0].get_cells()), &world.get_persons()[0].get_cells());
    ASSERT_NE(&(copied_world.get_persons()[1].get_cells()), &world.get_persons()[1].get_cells());

    // Evolve the world and check that the copied world has not evolved
    copied_world.get_persons()[0].migrate_to(work, {0});
    copied_world.get_persons()[1].migrate_to(home, {0});
    ASSERT_NE(copied_world.get_persons()[0].get_location().get_type(),
              world.get_persons()[0].get_location().get_type());
    ASSERT_NE(copied_world.get_persons()[1].get_location().get_type(),
              world.get_persons()[1].get_location().get_type());
}
