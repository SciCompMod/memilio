/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "abm/parameters.h"
#include "abm/person.h"
#include "abm/model.h"
#include "abm/virus_variant.h"
#include "abm_helpers.h"
#include "memilio/utils/parameter_distribution_wrapper.h"
#include "memilio/utils/parameter_distributions.h"
#include "random_number_test.h"

using TestModel = RandomNumberTest;

/**
 * @brief Test initialization of Model class.
 */
TEST_F(TestModel, init)
{
    auto model = mio::abm::Model(num_age_groups);

    // Verify that the model starts with exactly one location of type Cemetery.
    EXPECT_EQ(model.get_locations().size(), 1);
    EXPECT_EQ(model.get_locations()[0].get_type(), mio::abm::LocationType::Cemetery);

    // Verify that no persons are initialized in the model.
    EXPECT_THAT(model.get_persons(), testing::ElementsAre());
}

/**
 * @brief Test adding locations to the Model class.
 */
TEST_F(TestModel, addLocation)
{
    auto model      = mio::abm::Model(num_age_groups);
    auto school_id1 = model.add_location(mio::abm::LocationType::School);
    auto school_id2 = model.add_location(mio::abm::LocationType::School);
    auto work_id    = model.add_location(mio::abm::LocationType::Work);
    auto home_id    = model.add_location(mio::abm::LocationType::Home);

    // Verify the unique IDs of added locations.
    EXPECT_EQ(school_id1.get(), 1u);
    EXPECT_EQ(school_id2.get(), 2u);

    // Retrieve added locations by their IDs and verify their types.
    auto& school1 = model.get_location(school_id1);
    auto& school2 = model.get_location(school_id2);
    auto& work    = model.get_location(work_id);
    auto& home    = model.get_location(home_id);

    // Count the number of School-type locations and verify.
    size_t count_schools = 0;
    for (auto& loc : model.get_locations()) {
        if (loc.get_type() == mio::abm::LocationType::School) {
            count_schools++;
        }
    }
    EXPECT_EQ(count_schools, 2);

    // Verify the location order within the model's internal list of locations.
    EXPECT_EQ(model.get_locations()[1], school1);
    EXPECT_EQ(model.get_locations()[2], school2);
    EXPECT_EQ(model.get_locations()[3], work);
    EXPECT_EQ(model.get_locations()[4], home);
}

/**
 * @brief Test adding persons to a specific location in the Model class.
 */
TEST_F(TestModel, addPerson)
{
    auto model    = mio::abm::Model(num_age_groups);
    auto location = model.add_location(mio::abm::LocationType::School);

    model.add_person(location, age_group_15_to_34);
    model.add_person(location, age_group_35_to_59);

    // Verify the number of persons in the model and their respective age groups.
    EXPECT_EQ(model.get_persons().size(), 2);
    EXPECT_EQ(model.get_person(0).get_age(), age_group_15_to_34);
    EXPECT_EQ(model.get_person(1).get_age(), age_group_35_to_59);
    // Verify the number of persons in the model and their respective age groups.
    EXPECT_EQ(model.get_persons().size(), 2);
    EXPECT_EQ(model.get_person(0).get_age(), age_group_15_to_34);
    EXPECT_EQ(model.get_person(1).get_age(), age_group_35_to_59);
}

/**
 * @brief Test combined subpopulation count by location type in the Model class.
 */
TEST_F(TestModel, getSubpopulationCombined)
{
    auto t       = mio::abm::TimePoint(0);
    auto model   = mio::abm::Model(num_age_groups);
    auto school1 = model.add_location(mio::abm::LocationType::School);
    auto school2 = model.add_location(mio::abm::LocationType::School);
    auto school3 = model.add_location(mio::abm::LocationType::School);
    auto home1   = model.add_location(mio::abm::LocationType::Home);

    // Add persons to these locations with various infection states.
    add_test_person(model, school1, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);
    add_test_person(model, school1, age_group_15_to_34, mio::abm::InfectionState::Susceptible);
    add_test_person(model, school2, age_group_15_to_34, mio::abm::InfectionState::Susceptible);
    add_test_person(model, school2, age_group_15_to_34, mio::abm::InfectionState::Susceptible);
    add_test_person(model, school3, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);
    add_test_person(model, home1, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);

    // Verify the count of susceptible persons across all School locations.
    EXPECT_EQ(model.get_subpopulation_combined_per_location_type(t, mio::abm::InfectionState::Susceptible,
                                                                 mio::abm::LocationType::School),
              3);
    // Verify the count of persons with no symptoms across all School locations.
    EXPECT_EQ(model.get_subpopulation_combined_per_location_type(t, mio::abm::InfectionState::InfectedNoSymptoms,
                                                                 mio::abm::LocationType::School),
              2);
    // Verify the total count of persons with no symptoms across all locations.
    EXPECT_EQ(model.get_subpopulation_combined(t, mio::abm::InfectionState::InfectedNoSymptoms), 3);
}

/**
 * @brief Test finding a location assigned to a person in the Model class.
 */
TEST_F(TestModel, findLocation)
{
    // Create a model and add different location types.
    auto model      = mio::abm::Model(num_age_groups);
    model.get_rng() = this->get_rng();

    auto home_id   = model.add_location(mio::abm::LocationType::Home);
    auto school_id = model.add_location(mio::abm::LocationType::School);
    auto work_id   = model.add_location(mio::abm::LocationType::Work);

    // Add a person to the model and assign them to multiple locations.
    auto person_id = add_test_person(model, home_id);
    auto& person   = model.get_person(person_id);
    person.set_assigned_location(mio::abm::LocationType::Home, home_id);
    person.set_assigned_location(mio::abm::LocationType::Work, work_id);
    person.set_assigned_location(mio::abm::LocationType::School, school_id);

    // Verify that the find_location method correctly identifies each assigned location.
    EXPECT_EQ(model.find_location(mio::abm::LocationType::Work, person_id), work_id);
    EXPECT_EQ(model.find_location(mio::abm::LocationType::School, person_id), school_id);
    EXPECT_EQ(model.find_location(mio::abm::LocationType::Home, person_id), home_id);

    // Check that the method also works with a constant reference to the model.
    auto&& model_test = std::as_const(model);
    EXPECT_EQ(model_test.find_location(mio::abm::LocationType::Work, person_id), work_id);
    EXPECT_EQ(model_test.find_location(mio::abm::LocationType::School, person_id), school_id);
    EXPECT_EQ(model_test.find_location(mio::abm::LocationType::Home, person_id), home_id);
}

/**
 * @brief Test state transitions during a time step in the Model class.
 */
TEST_F(TestModel, evolveStateTransition)
{
    using testing::Return;

    auto t          = mio::abm::TimePoint(0);
    auto dt         = mio::abm::hours(1);
    auto model      = mio::abm::Model(num_age_groups);
    model.get_rng() = this->get_rng();

    // Setup incubation and infection period parameters to prevent state transitions within one hour. p1 and p3 don't transition.
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::LogNormalDistribution<double>>>> mock_logNorm_dist;
    EXPECT_CALL(mock_logNorm_dist.get_mock(), invoke).WillRepeatedly(testing::Return(2 * dt.days()));

    // Add locations and persons to the model with different initial infection states.
    auto location1 = model.add_location(mio::abm::LocationType::School);
    auto location2 = model.add_location(mio::abm::LocationType::Work);
    add_test_person(model, location1, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);
    add_test_person(model, location1, age_group_15_to_34, mio::abm::InfectionState::Susceptible);
    add_test_person(model, location2, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms);

    auto& p1 = model.get_persons()[0];
    auto& p2 = model.get_persons()[1];
    auto& p3 = model.get_persons()[2];

    // Assign persons to their respective locations.
    p1.set_assigned_location(mio::abm::LocationType::School, location1);
    p2.set_assigned_location(mio::abm::LocationType::School, location1);
    p3.set_assigned_location(mio::abm::LocationType::Work, location2);

    // Setup mock so p2 becomes infected
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.0));

    model.evolve(t, dt);

    // Verify the state transitions.
    EXPECT_EQ(p1.get_infection_state(t + dt), mio::abm::InfectionState::InfectedNoSymptoms);
    EXPECT_EQ(p2.get_infection_state(t + dt), mio::abm::InfectionState::Exposed);
    EXPECT_EQ(p3.get_infection_state(t + dt), mio::abm::InfectionState::InfectedSymptoms);
}

/**
 * @brief Test mobility rule-based transitions during a time step in the Model class.
 */
TEST_F(TestModel, evolveMobilityRules)
{
    using testing::Return;

    auto t          = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto dt         = mio::abm::hours(1);
    auto model      = mio::abm::Model(num_age_groups);
    model.get_rng() = this->get_rng();

    // Setup infection period parameters to prevent state transitions within one hour. p1 doesn't transition.
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::LogNormalDistribution<double>>>> mock_logNorm_dist;
    EXPECT_CALL(mock_logNorm_dist.get_mock(), invoke).WillRepeatedly(testing::Return(2 * dt.days()));
    model.parameters.get<mio::abm::AgeGroupGotoSchool>().set_multiple({age_group_5_to_14}, true);
    model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({age_group_15_to_34, age_group_35_to_59}, true);

    auto home_id   = model.add_location(mio::abm::LocationType::Home);
    auto school_id = model.add_location(mio::abm::LocationType::School);
    auto work_id   = model.add_location(mio::abm::LocationType::Work);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
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

    auto pid2 = add_test_person(model, home_id, age_group_5_to_14, mio::abm::InfectionState::Susceptible, t);
    auto pid1 = add_test_person(model, home_id, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms, t);

    auto& p1 = model.get_person(pid1);
    auto& p2 = model.get_person(pid2);

    p1.set_assigned_location(mio::abm::LocationType::School, school_id);
    p2.set_assigned_location(mio::abm::LocationType::School, school_id);
    p1.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p2.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p1.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p2.set_assigned_location(mio::abm::LocationType::Home, home_id);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.)); //no state transitions

    // Evolve the model over one time step and verify the location transitions.
    model.evolve(t, dt);
    EXPECT_EQ(p1.get_location(), work_id);
    EXPECT_EQ(p2.get_location(), school_id);

    // Verify the number of persons at each location.
    EXPECT_EQ(model.get_number_persons(school_id), 1);
    EXPECT_EQ(model.get_number_persons(work_id), 1);
}

/**
 * @brief Test the evolution of mobility trips within the Model class.
 */
TEST_F(TestModel, evolveMobilityTrips)
{
    using testing::Return;

    auto t     = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto dt    = mio::abm::hours(2);
    auto model = mio::abm::Model(num_age_groups);
    mio::ParameterDistributionConstant constant(2 * dt.days());
    //setup so p1-p5 don't do transition
    model.parameters
        .get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        mio::ParameterDistributionWrapper(constant);
    model.parameters
        .get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        mio::ParameterDistributionWrapper(constant);
    model.parameters
        .get<mio::abm::TimeInfectedSevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        mio::ParameterDistributionWrapper(constant);
    model.parameters
        .get<mio::abm::TimeInfectedSevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        mio::ParameterDistributionWrapper(constant);

    // Add different location types to the model.
    auto home_id     = model.add_location(mio::abm::LocationType::Home);
    auto event_id    = model.add_location(mio::abm::LocationType::SocialEvent);
    auto work_id     = model.add_location(mio::abm::LocationType::Work);
    auto hospital_id = model.add_location(mio::abm::LocationType::Hospital);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
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
        .WillRepeatedly(testing::Return(0.8)); // this forces p1 and p3 to recover

    // Create persons with various infection states and assign them to multiple locations.
    auto pid1 = add_test_person(model, home_id, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms, t);
    auto pid2 = add_test_person(model, home_id, age_group_15_to_34, mio::abm::InfectionState::Susceptible, t);
    auto pid3 = add_test_person(model, home_id, age_group_15_to_34, mio::abm::InfectionState::InfectedSevere, t);
    auto pid4 = add_test_person(model, hospital_id, age_group_15_to_34, mio::abm::InfectionState::Recovered, t);
    auto pid5 = add_test_person(model, home_id, age_group_15_to_34, mio::abm::InfectionState::Susceptible, t);

    // Assign persons to locations for trips.
    auto& p1 = model.get_person(pid1);
    auto& p2 = model.get_person(pid2);
    auto& p3 = model.get_person(pid3);
    auto& p4 = model.get_person(pid4);
    auto& p5 = model.get_person(pid5);

    p1.set_assigned_location(mio::abm::LocationType::SocialEvent, event_id);
    p2.set_assigned_location(mio::abm::LocationType::SocialEvent, event_id);
    p1.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p2.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p1.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p2.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p3.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p4.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p3.set_assigned_location(mio::abm::LocationType::Hospital, hospital_id);
    p4.set_assigned_location(mio::abm::LocationType::Hospital, hospital_id);
    p5.set_assigned_location(mio::abm::LocationType::SocialEvent, event_id);
    p5.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p5.set_assigned_location(mio::abm::LocationType::Home, home_id);

    // Set trips for persons between assigned locations.
    mio::abm::TripList& data = model.get_trip_list();
    mio::abm::Trip trip1(p1.get_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), work_id, home_id);
    mio::abm::Trip trip2(p2.get_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), event_id, home_id);
    mio::abm::Trip trip3(p5.get_id(), mio::abm::TimePoint(0) + mio::abm::hours(9), event_id, home_id);
    data.add_trip(trip1);
    data.add_trip(trip2);
    data.add_trip(trip3);

    // Set trips to use weekday trips on weekends.
    data.use_weekday_trips_on_weekend();
    // Set trips to use weekday trips on weekends.
    data.use_weekday_trips_on_weekend();

    // Mock the distribution to prevent infectionsin the test.
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.));

    model.evolve(t, dt);

    // Verify all the mobility transitions are correct.
    EXPECT_EQ(p1.get_location(), work_id);
    EXPECT_EQ(p2.get_location(), event_id);
    EXPECT_EQ(p3.get_location(), hospital_id);
    EXPECT_EQ(p4.get_location(), home_id);
    EXPECT_EQ(p5.get_location(), event_id);
    EXPECT_EQ(model.get_number_persons(event_id), 2);
    EXPECT_EQ(model.get_number_persons(work_id), 1);
    EXPECT_EQ(model.get_number_persons(home_id), 1);
    EXPECT_EQ(model.get_number_persons(hospital_id), 1);

    // Move all persons back to their home location to prepare for weekend trips.
    model.change_location(p1.get_id(), home_id);
    model.change_location(p3.get_id(), home_id);
    model.change_location(p2.get_id(), home_id);
    model.change_location(p5.get_id(), home_id);

    // Update the time to the weekend and reset the trip index.
    t = mio::abm::TimePoint(0) + mio::abm::days(6) + mio::abm::hours(8);
    model.get_trip_list().reset_index();

    // Evolve the model again to verify the weekend behavior.
    model.evolve(t, dt);

    EXPECT_EQ(p1.get_location(), work_id);
    EXPECT_EQ(p2.get_location(), event_id);
    EXPECT_EQ(p3.get_location(), home_id);
    EXPECT_EQ(p4.get_location(), home_id);
    EXPECT_EQ(p5.get_location(), event_id);
    EXPECT_EQ(model.get_number_persons(event_id), 2);
    EXPECT_EQ(model.get_number_persons(work_id), 1);
    EXPECT_EQ(model.get_number_persons(home_id), 2);

    // Add additional weekend trips for further verification.
    bool weekend = true;
    mio::abm::Trip tripweekend1(p1.get_id(), mio::abm::TimePoint(0) + mio::abm::days(6) + mio::abm::hours(10),
                                event_id);
    mio::abm::Trip tripweekend2(p2.get_id(), mio::abm::TimePoint(0) + mio::abm::days(6) + mio::abm::hours(10), home_id);
    mio::abm::Trip tripweekend3(p5.get_id(), mio::abm::TimePoint(0) + mio::abm::days(6) + mio::abm::hours(10), work_id);
    data.add_trip(tripweekend1, weekend);
    data.add_trip(tripweekend2, weekend);
    data.add_trip(tripweekend3, weekend);

    // Advance time and evolve the model to check location transitions during the weekend.
    t += mio::abm::hours(1);

    model.evolve(t, dt);

    EXPECT_EQ(p1.get_location(), event_id);
    EXPECT_EQ(p2.get_location(), home_id);
    EXPECT_EQ(p3.get_location(), home_id);
    EXPECT_EQ(p4.get_location(), home_id);
    EXPECT_EQ(p5.get_location(), work_id);
    EXPECT_EQ(model.get_number_persons(event_id), 1);
    EXPECT_EQ(model.get_number_persons(work_id), 1);
    EXPECT_EQ(model.get_number_persons(home_id), 3);
}

#ifndef MEMILIO_ENABLE_OPENMP // TODO: Test can fail with parallel execution of mobility, as the capacity is not taken into account correctly at the moment (c. f. issue #640)
/**
 * @brief Test that a location correctly enforces its capacity constraint.
 */
TEST_F(TestModel, reachCapacity)
{
    using testing::Return;

    // Initialize time and model.
    auto t          = mio::abm::TimePoint{mio::abm::hours(8).seconds()};
    auto dt         = mio::abm::hours(1);
    auto model      = mio::abm::Model(num_age_groups);
    model.get_rng() = this->get_rng();
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;

    auto home_id   = model.add_location(mio::abm::LocationType::Home);
    auto school_id = model.add_location(mio::abm::LocationType::School);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
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

    // Create two persons with different infection states.
    auto p1 = add_test_person(model, home_id, age_group_5_to_14);
    auto p2 = add_test_person(model, home_id, age_group_5_to_14);

    // Assign both persons to School and Home.
    model.get_person(p1).set_assigned_location(mio::abm::LocationType::School, school_id);
    model.get_person(p2).set_assigned_location(mio::abm::LocationType::School, school_id);
    model.get_person(p1).set_assigned_location(mio::abm::LocationType::Home, home_id);
    model.get_person(p2).set_assigned_location(mio::abm::LocationType::Home, home_id);

    // Set the capacity of the school to 1 person with a distance requirement of 66.
    model.get_location(school_id).set_capacity(1, 66);

    model.evolve(t, dt);

    // Verify that only one person is at the school, while the other remains at home due to capacity constraints.
    EXPECT_EQ(model.get_person(p1).get_location(), school_id);
    EXPECT_EQ(model.get_person(p2).get_location(), home_id); // p2 should not be able to enter the school
    EXPECT_EQ(model.get_number_persons(school_id), 1);
    EXPECT_EQ(model.get_number_persons(home_id), 1);
}
#endif

/**
 * @brief Test that dead persons remain in the cemetery and can't be moved by scheduled trips.
 */
TEST_F(TestModel, checkMobilityOfDeadPerson)
{
    using testing::Return;
    auto t     = mio::abm::TimePoint(0);
    auto dt    = mio::abm::days(1);
    auto model = mio::abm::Model(num_age_groups);
    model.parameters
        .get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 1.;
    model.parameters
        .get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 1.;
    // Time to go from severe to critical infection is 1 day (dt).
    mio::ParameterDistributionConstant constant1(dt.days());
    model.parameters
        .get<mio::abm::TimeInfectedSevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] =
        mio::ParameterDistributionWrapper(constant1);
    // Time to go from critical infection to dead state is 1/2 day (0.5 * dt).
    mio::ParameterDistributionConstant constant2(0.5 * dt.days());
    model.parameters
        .get<mio::abm::TimeInfectedCriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] =
        mio::ParameterDistributionWrapper(constant2);

    auto home_id     = model.add_location(mio::abm::LocationType::Home);
    auto work_id     = model.add_location(mio::abm::LocationType::Work);
    auto icu_id      = model.add_location(mio::abm::LocationType::ICU);
    auto hospital_id = model.add_location(mio::abm::LocationType::Hospital);
    // Create a person that is dead at time t
    add_test_person(model, icu_id, age_group_60_to_79, mio::abm::InfectionState::Dead, t);
    // Create a person that is severe at hospital and will be dead at time t + dt
    add_test_person(model, hospital_id, age_group_60_to_79, mio::abm::InfectionState::Dead, t + dt);

    auto& p_dead   = model.get_persons()[0];
    auto& p_severe = model.get_persons()[1];

    p_dead.set_assigned_location(mio::abm::LocationType::ICU, icu_id);
    p_dead.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p_dead.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p_severe.set_assigned_location(mio::abm::LocationType::Hospital, hospital_id);
    p_severe.set_assigned_location(mio::abm::LocationType::ICU, icu_id);
    p_severe.set_assigned_location(mio::abm::LocationType::Home, home_id);

    // Add trip to see if a dead person can change location outside of cemetery by scheduled trips
    mio::abm::TripList& trip_list = model.get_trip_list();
    mio::abm::Trip trip1(p_dead.get_id(), mio::abm::TimePoint(0) + mio::abm::hours(2), work_id, home_id);
    mio::abm::Trip trip2(p_dead.get_id(), mio::abm::TimePoint(0) + mio::abm::hours(3), home_id, icu_id);
    mio::abm::Trip trip3(p_severe.get_id(), mio::abm::TimePoint(0) + mio::abm::hours(3), home_id, icu_id);
    trip_list.add_trip(trip1);
    trip_list.add_trip(trip2);
    trip_list.add_trip(trip3);

    // Check the dead person got burried and the severely infected person starts in Hospital
    model.evolve(t, dt);
    EXPECT_EQ(model.get_location(p_dead.get_id()).get_type(), mio::abm::LocationType::Cemetery);
    EXPECT_EQ(p_severe.get_infection_state(t), mio::abm::InfectionState::InfectedSevere);
    EXPECT_EQ(model.get_location(p_severe.get_id()).get_type(), mio::abm::LocationType::Hospital);

    // Check the dead person is still in Cemetery and the severely infected person dies and got burried
    model.evolve(t + dt, dt);
    EXPECT_EQ(model.get_location(p_dead.get_id()).get_type(), mio::abm::LocationType::Cemetery);
    EXPECT_EQ(p_severe.get_infection_state(t + dt), mio::abm::InfectionState::Dead);
    EXPECT_EQ(model.get_location(p_severe.get_id()).get_type(), mio::abm::LocationType::Cemetery);
}

using TestModelTestingCriteria = RandomNumberTest;

/**
 * @brief Test the adding, updating, and running of testing schemes within the model.
 */
TEST_F(TestModelTestingCriteria, testAddingAndUpdatingAndRunningTestingSchemes)
{
    auto model = mio::abm::Model(num_age_groups);
    // make sure the infected person stay in Infected long enough
    mio::ParameterDistributionConstant constant(100.);
    model.parameters.get<mio::abm::TimeInfectedSymptomsToRecovered>()[{mio::abm::VirusVariant(0), age_group_15_to_34}] =
        mio::ParameterDistributionWrapper(constant);
    model.parameters.get<mio::abm::TimeInfectedSymptomsToSevere>()[{mio::abm::VirusVariant(0), age_group_15_to_34}] =
        mio::ParameterDistributionWrapper(constant);

    auto home_id = model.add_location(mio::abm::LocationType::Home);
    auto work_id = model.add_location(mio::abm::LocationType::Work);
    auto& work   = model.get_location(work_id);

    auto current_time = mio::abm::TimePoint(0);

    auto test_time = mio::abm::minutes(30);
    // Add a person to the model with an infection state that requires testing.
    // Since tests are performed before current_time, the InfectionState of the Person has to take into account test_time
    auto pid        = add_test_person(model, home_id, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms,
                                      current_time - test_time);
    auto& person    = model.get_person(pid);
    auto rng_person = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
    person.set_assigned_location(mio::abm::LocationType::Home, home_id);
    person.set_assigned_location(mio::abm::LocationType::Work, work_id);

    auto validity_period       = mio::abm::days(1);
    const auto start_date      = mio::abm::TimePoint(20);
    const auto end_date        = mio::abm::TimePoint(60 * 60 * 24 * 3);
    const auto probability     = 1.0;
    const auto test_params_pcr = mio::abm::TestParameters{0.9, 0.99, test_time, mio::abm::TestType::Generic};

    auto testing_criteria = mio::abm::TestingCriteria();
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedSymptoms);
    testing_criteria.add_infection_state(mio::abm::InfectionState::InfectedNoSymptoms);

    auto testing_scheme =
        mio::abm::TestingScheme(testing_criteria, validity_period, start_date, end_date, test_params_pcr, probability);

    model.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme);
    EXPECT_EQ(model.get_testing_strategy().run_strategy(rng_person, person, work, current_time),
              true); // no active testing scheme -> person can enter
    current_time = mio::abm::TimePoint(30);
    model.get_testing_strategy().update_activity_status(current_time);
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::Exactly(5))
        .WillOnce(testing::Return(0.7)) // Person complies with testing
        .WillOnce(testing::Return(0.5)) // Probability for testing (is performed)
        .WillOnce(testing::Return(0.4)) // Test result is positive
        .WillOnce(testing::Return(0.0)) // Draw for isolation compliance (doesn't matter in this test)
        .WillOnce(
            testing::Return(0.7)); // Person complies with testing (even though there is not testing strategy left)
    EXPECT_EQ(model.get_testing_strategy().run_strategy(rng_person, person, work, current_time),
              false); // Testing scheme active and restricts entry

    // Try to re-add the same testing scheme and confirm it doesn't duplicate, then remove it.
    model.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work,
                                                    testing_scheme); //doesn't get added because of == operator
    model.get_testing_strategy().remove_testing_scheme(mio::abm::LocationType::Work, testing_scheme);
    EXPECT_EQ(model.get_testing_strategy().run_strategy(rng_person, person, work, current_time),
              true); // no more testing_schemes
}

/**
 * @brief Test to validate the parameter constraints within the model.
 */
TEST_F(TestModel, checkParameterConstraints)
{
    mio::set_log_level(mio::LogLevel::critical); // Suppress logging of errors since they are expected here
    auto model  = mio::abm::Model(num_age_groups);
    auto params = model.parameters;

    mio::ParameterDistributionLogNormal log_norm1(1., 0.5);
    mio::ParameterDistributionLogNormal log_norm2(2., 0.5);
    mio::ParameterDistributionLogNormal log_norm3(3., 0.5);
    mio::ParameterDistributionLogNormal log_norm4(4., 0.5);
    mio::ParameterDistributionLogNormal log_norm5(5., 0.5);
    mio::ParameterDistributionLogNormal log_norm6(6., 0.5);
    mio::ParameterDistributionLogNormal log_norm7(7., 0.5);
    mio::ParameterDistributionLogNormal log_norm8(8., 0.5);
    mio::ParameterDistributionLogNormal log_norm9(9., 0.5);
    params.get<mio::abm::IncubationPeriod>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm1);
    params.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm2);
    params.get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm3);
    params.get<mio::abm::TimeInfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm4);
    params.get<mio::abm::TimeInfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm5);
    params.get<mio::abm::TimeInfectedSevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm6);
    params.get<mio::abm::TimeInfectedSevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm7);
    params.get<mio::abm::TimeInfectedCriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm8);
    params.get<mio::abm::TimeInfectedCriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm9);
    params.get<mio::abm::DetectInfection>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = 0.3;
    params.get<mio::abm::GotoWorkTimeMinimum>()[age_group_35_to_59]                               = mio::abm::hours(4);
    params.get<mio::abm::GotoWorkTimeMaximum>()[age_group_35_to_59]                               = mio::abm::hours(8);
    params.get<mio::abm::GotoSchoolTimeMinimum>()[age_group_0_to_4]                               = mio::abm::hours(3);
    params.get<mio::abm::GotoSchoolTimeMaximum>()[age_group_0_to_4]                               = mio::abm::hours(6);
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::Community]                         = 0.5;
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::FFP2]                              = 0.6;
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::Surgical]                          = 0.7;
    params.get<mio::abm::LockdownDate>() = mio::abm::TimePoint(0);
    ASSERT_EQ(params.check_constraints(), false);

    mio::ParameterDistributionLogNormal log_normm1(-1., 0.5);
    mio::ParameterDistributionLogNormal log_normm2(-2., 0.5);
    mio::ParameterDistributionLogNormal log_normm3(-3., 0.5);
    mio::ParameterDistributionLogNormal log_normm4(-4., 0.5);
    mio::ParameterDistributionLogNormal log_normm5(-5., 0.5);
    mio::ParameterDistributionLogNormal log_normm6(-6., 0.5);
    mio::ParameterDistributionLogNormal log_normm7(-7., 0.5);
    mio::ParameterDistributionLogNormal log_normm8(-8., 0.5);
    mio::ParameterDistributionLogNormal log_normm9(-9., 0.5);

    params.get<mio::abm::IncubationPeriod>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_normm1);
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::IncubationPeriod>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm1);
    params.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_normm2);
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm2);
    params.get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_normm3);
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm3);
    params.get<mio::abm::TimeInfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_normm4);
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::TimeInfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm4);
    params.get<mio::abm::TimeInfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_normm5);
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::TimeInfectedSymptomsToSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm5);
    params.get<mio::abm::TimeInfectedSevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_normm6);
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::TimeInfectedSevereToCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm6);
    params.get<mio::abm::TimeInfectedSevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_normm7);
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::TimeInfectedSevereToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm7);
    params.get<mio::abm::TimeInfectedCriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_normm8);
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::TimeInfectedCriticalToDead>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm8);
    params.get<mio::abm::TimeInfectedCriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_normm9);
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::TimeInfectedCriticalToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] =
        mio::ParameterDistributionWrapper(log_norm9);
    params.get<mio::abm::DetectInfection>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = 1.1;
    ASSERT_EQ(params.check_constraints(), true);
    params.get<mio::abm::DetectInfection>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}] = 0.3;

    params.get<mio::abm::GotoWorkTimeMinimum>()[age_group_35_to_59] = mio::abm::hours(30);
    EXPECT_TRUE(params.check_constraints());
    params.get<mio::abm::GotoWorkTimeMinimum>()[age_group_35_to_59] = mio::abm::hours(4);
    params.get<mio::abm::GotoWorkTimeMaximum>()[age_group_35_to_59] = mio::abm::hours(30);
    EXPECT_TRUE(params.check_constraints());
    params.get<mio::abm::GotoWorkTimeMaximum>()[age_group_35_to_59] = mio::abm::hours(8);
    params.get<mio::abm::GotoSchoolTimeMinimum>()[age_group_0_to_4] = mio::abm::hours(30);
    EXPECT_TRUE(params.check_constraints());
    params.get<mio::abm::GotoSchoolTimeMinimum>()[age_group_0_to_4] = mio::abm::hours(3);
    params.get<mio::abm::GotoSchoolTimeMaximum>()[age_group_0_to_4] = mio::abm::hours(30);
    EXPECT_TRUE(params.check_constraints());
    params.get<mio::abm::GotoSchoolTimeMaximum>()[age_group_0_to_4] = mio::abm::hours(6);

    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::Community] = 1.2;
    EXPECT_TRUE(params.check_constraints());
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::Community] = 0.5;
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::FFP2]      = 1.2;
    EXPECT_TRUE(params.check_constraints());
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::FFP2]     = 0.6;
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::Surgical] = 1.2;
    EXPECT_TRUE(params.check_constraints());
    params.get<mio::abm::MaskProtection>()[mio::abm::MaskType::Surgical] = 0.7;

    params.get<mio::abm::LockdownDate>() = mio::abm::TimePoint(-2);
    EXPECT_TRUE(params.check_constraints());
    mio::set_log_level(mio::LogLevel::warn);
}

/**
 * @brief Test the enforcement of NPIs (Non-Pharmaceutical Interventions) on mobility rules.
 */
TEST_F(TestModel, mobilityRulesWithAppliedNPIs)
{
    using testing::Return;
    // Test when the NPIs are applied, people can enter targeted location if they comply to the rules.
    auto t          = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto dt         = mio::abm::hours(1);
    auto test_time  = mio::abm::minutes(30);
    auto model      = mio::abm::Model(num_age_groups);
    model.get_rng() = this->get_rng();

    mio::ParameterDistributionConstant constant(2 * dt.days());
    model.parameters
        .get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        mio::ParameterDistributionWrapper(constant);
    model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({age_group_15_to_34, age_group_35_to_59}, true);
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;

    auto home_id   = model.add_location(mio::abm::LocationType::Home);
    auto work_id   = model.add_location(mio::abm::LocationType::Work);
    auto school_id = model.add_location(mio::abm::LocationType::School);
    auto& work     = model.get_location(work_id);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(16))
        .WillOnce(testing::Return(0.8)) // draw random work group
        .WillOnce(testing::Return(0.8)) // draw random school group
        .WillOnce(testing::Return(0.8)) // draw random work hour
        .WillOnce(testing::Return(0.8)) // draw random school hour
        .WillOnce(testing::Return(0.8)) // draw random work group
        .WillOnce(testing::Return(0.8)) // draw random school group
        .WillOnce(testing::Return(0.8)) // draw random work hour
        .WillOnce(testing::Return(0.8)) // draw random school hour
        .WillOnce(testing::Return(0.8)) // draw random work group
        .WillOnce(testing::Return(0.8)) // draw random school group
        .WillOnce(testing::Return(0.8)) // draw random work hour
        .WillOnce(testing::Return(0.8)) // draw random school hour
        .WillOnce(testing::Return(0.8)) // draw random work group
        .WillOnce(testing::Return(0.8)) // draw random school group
        .WillOnce(testing::Return(0.8)) // draw random work hour
        .WillOnce(testing::Return(0.8)) // draw random school hour
        .WillRepeatedly(testing::Return(0.9)); // draw that satisfies all pre-conditions of NPIs

    // Since tests are performed before t, the InfectionState of all the Person have to take into account test_time
    auto p_id_compliant_go_to_work =
        add_test_person(model, home_id, age_group_15_to_34, mio::abm::InfectionState::Susceptible, t - test_time);
    auto p_id_compliant_go_to_school =
        add_test_person(model, home_id, age_group_5_to_14, mio::abm::InfectionState::Susceptible, t - test_time);
    auto p_id_no_mask =
        add_test_person(model, home_id, age_group_15_to_34, mio::abm::InfectionState::Susceptible, t - test_time);
    auto p_id_no_test      = add_test_person(model, home_id, age_group_15_to_34,
                                             mio::abm::InfectionState::InfectedNoSymptoms, t - test_time);
    auto p_id_no_isolation = add_test_person(model, home_id, age_group_15_to_34,
                                             mio::abm::InfectionState::InfectedNoSymptoms, t - test_time);

    auto& p_compliant_go_to_work   = model.get_person(p_id_compliant_go_to_work);
    auto& p_compliant_go_to_school = model.get_person(p_id_compliant_go_to_school);
    auto& p_no_mask                = model.get_person(p_id_no_mask);
    auto& p_no_test                = model.get_person(p_id_no_test);
    auto& p_no_isolation           = model.get_person(p_id_no_isolation);

    p_compliant_go_to_work.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p_compliant_go_to_work.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p_compliant_go_to_work.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p_compliant_go_to_school.set_assigned_location(mio::abm::LocationType::School, school_id);
    p_compliant_go_to_school.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p_no_mask.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p_no_mask.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p_no_test.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p_no_test.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p_no_isolation.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p_no_isolation.set_assigned_location(mio::abm::LocationType::Home, home_id);

    auto testing_criteria = mio::abm::TestingCriteria(
        {}, {mio::abm::InfectionState::InfectedSymptoms, mio::abm::InfectionState::InfectedNoSymptoms});
    const auto start_date        = mio::abm::TimePoint(0);
    const auto end_date          = mio::abm::TimePoint(60 * 60 * 24 * 3);
    const auto probability       = 1;
    const auto test_params       = mio::abm::TestParameters{0.99, 0.99, test_time, mio::abm::TestType::Generic};
    const auto testing_frequency = mio::abm::days(1);

    auto testing_scheme =
        mio::abm::TestingScheme(testing_criteria, testing_frequency, start_date, end_date, test_params, probability);
    model.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.));

    work.set_required_mask(mio::abm::MaskType::FFP2);
    p_no_mask.set_compliance(mio::abm::InterventionType::Mask, 0.4);
    p_no_test.set_compliance(mio::abm::InterventionType::Testing, 0.4);
    p_no_isolation.set_compliance(mio::abm::InterventionType::Isolation, 0.4);

    model.evolve(t, dt);

    // The complied person is allowed to be at work and wear the required mask
    EXPECT_EQ(p_compliant_go_to_work.get_location(), work_id);
    EXPECT_EQ(p_compliant_go_to_work.get_mask().get_type(), mio::abm::MaskType::FFP2);

    // The complied person is allowed to be at school and don't wear mask
    EXPECT_EQ(p_compliant_go_to_school.get_location(), school_id);
    EXPECT_EQ(p_compliant_go_to_school.get_mask().get_type(), mio::abm::MaskType::None);

    // The person, who does not wear mask, is not allowed to be in location
    EXPECT_EQ(p_no_mask.get_mask().get_type(), mio::abm::MaskType::None);
    EXPECT_NE(p_no_mask.get_location(), work_id);

    // The person, who does not want test, is not allowed to be in location
    EXPECT_NE(p_no_test.get_location(), work_id);

    // The person does not want to isolate when the test is positive
    EXPECT_FALSE(p_no_isolation.is_in_quarantine(t, model.parameters));
}

/**
 * @brief Test the enforcement of NPIs (Non-Pharmaceutical Interventions) on trips.
 */
TEST_F(TestModel, mobilityTripWithAppliedNPIs)
{
    using testing::Return;
    // Test when the NPIs are applied, people can enter targeted location if they comply to the rules.
    auto t          = mio::abm::TimePoint(0) + mio::abm::hours(8);
    auto dt         = mio::abm::hours(1);
    auto test_time  = mio::abm::minutes(30);
    auto model      = mio::abm::Model(num_age_groups);
    model.get_rng() = this->get_rng();

    mio::ParameterDistributionConstant constant(2 * dt.days());
    model.parameters
        .get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        mio::ParameterDistributionWrapper(constant);
    model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({age_group_15_to_34, age_group_35_to_59}, true);
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;

    auto home_id   = model.add_location(mio::abm::LocationType::Home);
    auto work_id   = model.add_location(mio::abm::LocationType::Work);
    auto school_id = model.add_location(mio::abm::LocationType::School);
    auto& work     = model.get_location(work_id);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(16))
        .WillOnce(testing::Return(0.8)) // draw random work group
        .WillOnce(testing::Return(0.8)) // draw random school group
        .WillOnce(testing::Return(0.8)) // draw random work hour
        .WillOnce(testing::Return(0.8)) // draw random school hour
        .WillOnce(testing::Return(0.8)) // draw random work group
        .WillOnce(testing::Return(0.8)) // draw random school group
        .WillOnce(testing::Return(0.8)) // draw random work hour
        .WillOnce(testing::Return(0.8)) // draw random school hour
        .WillOnce(testing::Return(0.8)) // draw random work group
        .WillOnce(testing::Return(0.8)) // draw random school group
        .WillOnce(testing::Return(0.8)) // draw random work hour
        .WillOnce(testing::Return(0.8)) // draw random school hour
        .WillOnce(testing::Return(0.8)) // draw random work group
        .WillOnce(testing::Return(0.8)) // draw random school group
        .WillOnce(testing::Return(0.8)) // draw random work hour
        .WillOnce(testing::Return(0.8)) // draw random school hour
        .WillRepeatedly(testing::Return(0.9)); // draw that satisfies all pre-conditions of NPIs

    // Since tests are performed before t, the InfectionState of all the Person have to take into account test_time
    auto p_id_compliant_go_to_work =
        add_test_person(model, home_id, age_group_15_to_34, mio::abm::InfectionState::Susceptible, t - test_time);
    auto p_id_compliant_go_to_school =
        add_test_person(model, home_id, age_group_5_to_14, mio::abm::InfectionState::Susceptible, t - test_time);
    auto p_id_no_mask =
        add_test_person(model, home_id, age_group_15_to_34, mio::abm::InfectionState::Susceptible, t - test_time);
    auto p_id_no_test      = add_test_person(model, home_id, age_group_15_to_34,
                                             mio::abm::InfectionState::InfectedNoSymptoms, t - test_time);
    auto p_id_no_isolation = add_test_person(model, home_id, age_group_15_to_34,
                                             mio::abm::InfectionState::InfectedNoSymptoms, t - test_time);

    auto& p_compliant_go_to_work   = model.get_person(p_id_compliant_go_to_work);
    auto& p_compliant_go_to_school = model.get_person(p_id_compliant_go_to_school);
    auto& p_no_mask                = model.get_person(p_id_no_mask);
    auto& p_no_test                = model.get_person(p_id_no_test);
    auto& p_no_isolation           = model.get_person(p_id_no_isolation);

    p_compliant_go_to_work.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p_compliant_go_to_work.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p_compliant_go_to_work.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p_compliant_go_to_school.set_assigned_location(mio::abm::LocationType::School, school_id);
    p_compliant_go_to_school.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p_no_mask.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p_no_mask.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p_no_test.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p_no_test.set_assigned_location(mio::abm::LocationType::Home, home_id);
    p_no_isolation.set_assigned_location(mio::abm::LocationType::Work, work_id);
    p_no_isolation.set_assigned_location(mio::abm::LocationType::Home, home_id);

    auto testing_criteria = mio::abm::TestingCriteria(
        {}, {mio::abm::InfectionState::InfectedSymptoms, mio::abm::InfectionState::InfectedNoSymptoms});
    const auto start_date        = mio::abm::TimePoint(0);
    const auto end_date          = mio::abm::TimePoint(60 * 60 * 24 * 3);
    const auto probability       = 1;
    const auto test_params       = mio::abm::TestParameters{0.99, 0.99, test_time, mio::abm::TestType::Generic};
    const auto testing_frequency = mio::abm::days(1);

    auto testing_scheme =
        mio::abm::TestingScheme(testing_criteria, testing_frequency, start_date, end_date, test_params, probability);
    model.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillRepeatedly(Return(1.));

    work.set_required_mask(mio::abm::MaskType::FFP2);
    p_no_mask.set_compliance(mio::abm::InterventionType::Mask, 0.4);
    p_no_test.set_compliance(mio::abm::InterventionType::Testing, 0.4);
    p_no_isolation.set_compliance(mio::abm::InterventionType::Isolation, 0.4);

    // Using trip list
    mio::abm::TripList& trip_list = model.get_trip_list();
    mio::abm::Trip trip1(p_compliant_go_to_work.get_id(), t, work_id, home_id);
    mio::abm::Trip trip2(p_compliant_go_to_school.get_id(), t, school_id, home_id);
    mio::abm::Trip trip3(p_no_mask.get_id(), t, work_id, home_id);
    mio::abm::Trip trip4(p_no_test.get_id(), t, work_id, home_id);
    mio::abm::Trip trip5(p_no_isolation.get_id(), t, work_id, home_id);
    trip_list.add_trip(trip1);
    trip_list.add_trip(trip2);
    trip_list.add_trip(trip3);
    trip_list.add_trip(trip4);
    trip_list.add_trip(trip5);
    model.use_mobility_rules(false);
    model.evolve(t, dt);

    // The complied person is allowed to be at work and wear the required mask
    EXPECT_EQ(p_compliant_go_to_work.get_location(), work_id);
    EXPECT_EQ(p_compliant_go_to_work.get_mask().get_type(), mio::abm::MaskType::FFP2);

    // The complied person is allowed to be at school and don't wear mask
    EXPECT_EQ(p_compliant_go_to_school.get_location(), school_id);
    EXPECT_EQ(p_compliant_go_to_school.get_mask().get_type(), mio::abm::MaskType::None);

    // The person, who does not wear mask, is not allowed to be in location
    EXPECT_EQ(p_no_mask.get_mask().get_type(), mio::abm::MaskType::None);
    EXPECT_NE(p_no_mask.get_location(), work_id);

    // The person, who does not want test, is not allowed to be in location
    EXPECT_NE(p_no_test.get_location(), work_id);

    // The person does not want to isolate when the test is positive
    EXPECT_FALSE(p_no_isolation.is_in_quarantine(t, model.parameters));
}
