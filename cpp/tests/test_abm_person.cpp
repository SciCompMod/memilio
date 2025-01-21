/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "abm/location_id.h"
#include "abm/model_functions.h"
#include "abm/location_type.h"
#include "abm/mobility_rules.h"
#include "abm/person.h"
#include "abm/time.h"
#include "abm_helpers.h"
#include "random_number_test.h"

#include <gtest/gtest.h>

using TestPerson = RandomNumberTest;

/**
 * @brief Test the initialization of a Person object with default properties.
 */
TEST_F(TestPerson, init)
{
    mio::abm::Location location(mio::abm::LocationType::Work, 7, num_age_groups);
    auto t      = mio::abm::TimePoint(0);
    auto person = mio::abm::Person(this->get_rng(), location.get_type(), location.get_id(), age_group_60_to_79);

    // Verify default state and location assignments.
    EXPECT_EQ(person.get_infection_state(t), mio::abm::InfectionState::Susceptible);
    EXPECT_EQ(person.get_location(), location.get_id());
    EXPECT_EQ(person.get_id(), mio::abm::PersonId::invalid_id());
}

/**
 * @brief Test that a Person's location can be changed correctly.
 */
TEST_F(TestPerson, change_location)
{
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location loc1(mio::abm::LocationType::PublicTransport, 1, 6, 1);
    mio::abm::Location loc2(mio::abm::LocationType::School, 2, num_age_groups);
    mio::abm::Location loc3(mio::abm::LocationType::PublicTransport, 3, 6, 2);
    auto person = make_test_person(this->get_rng(), home, age_group_0_to_4, mio::abm::InfectionState::Recovered);

    // Check that a person does not change location to its current location
    person.add_time_at_location(mio::abm::hours(1));
    EXPECT_FALSE(mio::abm::change_location(person, home));
    EXPECT_EQ(person.get_time_at_location(), mio::abm::hours(1));
    EXPECT_EQ(person.get_location(), home.get_id());

    // Change the location of the person a couple of times
    EXPECT_TRUE(mio::abm::change_location(person, loc1, mio::abm::TransportMode::Unknown, {0}));
    EXPECT_EQ(person.get_time_at_location(), mio::abm::TimeSpan(0));
    EXPECT_EQ(person.get_location(), loc1.get_id());
    EXPECT_EQ(person.get_last_transport_mode(), mio::abm::TransportMode::Unknown);

    EXPECT_TRUE(mio::abm::change_location(person, loc2, mio::abm::TransportMode::Walking, {0}));
    EXPECT_EQ(person.get_time_at_location(), mio::abm::TimeSpan(0));
    EXPECT_EQ(person.get_location(), loc2.get_id());
    EXPECT_EQ(person.get_last_transport_mode(), mio::abm::TransportMode::Walking);

    // Test changing location with cell indices.
    EXPECT_TRUE(mio::abm::change_location(person, loc3, mio::abm::TransportMode::Bike, {0, 1}));
    EXPECT_EQ(person.get_time_at_location(), mio::abm::TimeSpan(0));
    EXPECT_EQ(person.get_location(), loc3.get_id());
    EXPECT_EQ(person.get_last_transport_mode(), mio::abm::TransportMode::Bike);
    ASSERT_EQ(person.get_cells().size(), 2);
    EXPECT_EQ(person.get_cells()[0], 0u);
    EXPECT_EQ(person.get_cells()[1], 1u);
}

/**
 * @brief Test setting and retrieving assigned locations for a Person.
 */
TEST_F(TestPerson, setGetAssignedLocation)
{
    mio::abm::Location location(mio::abm::LocationType::Work, 2, num_age_groups);
    auto person = mio::abm::Person(this->get_rng(), location.get_type(), location.get_id(), age_group_35_to_59);
    // Assign and verify a location for the person.
    person.set_assigned_location(location.get_type(), location.get_id());
    EXPECT_EQ(person.get_assigned_location(mio::abm::LocationType::Work), mio::abm::LocationId(2));
    // Change the assigned location and verify.
    person.set_assigned_location(mio::abm::LocationType::Work, mio::abm::LocationId(4));
    EXPECT_EQ(person.get_assigned_location(mio::abm::LocationType::Work), mio::abm::LocationId(4));

    // Fuzzing: assign random valid LocationId values and verify correctness.
    for (int i = 0; i < 100; ++i) {
        auto random_id   = this->random_integer(0, std::numeric_limits<int>::max());
        auto random_type = this->random_integer(0, (int)mio::abm::LocationType::Count - 1);
        person.set_assigned_location((mio::abm::LocationType)random_type, mio::abm::LocationId(random_id));
        EXPECT_EQ(person.get_assigned_location((mio::abm::LocationType)random_type), mio::abm::LocationId(random_id));
    }

    // Boundary test cases: test with boundary LocationIds.
    person.set_assigned_location(mio::abm::LocationType::Work, mio::abm::LocationId(0));
    EXPECT_EQ(person.get_assigned_location(mio::abm::LocationType::Work), mio::abm::LocationId(0));

    person.set_assigned_location(mio::abm::LocationType::Work, mio::abm::LocationId(std::numeric_limits<int>::max()));
    EXPECT_EQ(person.get_assigned_location(mio::abm::LocationType::Work), mio::abm::LocationId(std::numeric_limits<int>::max()));
}

/**
 * @brief Test quarantine behavior and removal of isolation for a Person.
 */
TEST_F(TestPerson, quarantine)
{
    using testing::Return;
    auto test_params =
        mio::abm::TestParameters{1.01, 1.01, mio::abm::minutes(30), mio::abm::TestType::Generic}; //100% safe test

    auto infection_parameters = mio::abm::Parameters(num_age_groups);
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location work(mio::abm::LocationType::Work, 1, num_age_groups);

    // Setup rng mock so the person has a state transition to Recovered
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(4))
        .WillOnce(testing::Return(0.6)) // workgroup
        .WillOnce(testing::Return(0.6)) // schoolgroup
        .WillOnce(testing::Return(0.6)) // goto_work_hour
        .WillOnce(testing::Return(0.6)) // goto_school_hour
        .WillRepeatedly(testing::Return(1.0)); // ViralLoad draws

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::LogNormalDistribution<double>>>> mock_logNorm_dist;

    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(7);
    auto dt        = mio::abm::hours(1);
    EXPECT_CALL(mock_logNorm_dist.get_mock(), invoke)
        .Times(testing::AtLeast(1))
        .WillOnce(testing::Return(1.0)) // TimeInfectedNoSymptomsToSymptoms
        .WillOnce(testing::Return(1.0)) //IncubationTime
        .WillOnce(testing::Return(0.5 * dt.days())) // TimeInfectedSymptomsToRecovered
        .WillRepeatedly(testing::Return(1.0));

    infection_parameters.get<mio::abm::AgeGroupGotoSchool>().set_multiple({age_group_5_to_14}, true);
    infection_parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({age_group_15_to_34, age_group_35_to_59}, true);

    auto person     = make_test_person(this->get_rng(), home, age_group_35_to_59,
                                       mio::abm::InfectionState::InfectedSymptoms, t_morning, infection_parameters);
    auto rng_person = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), person);

    // Test quarantine when a person is tested and positive.
    person.get_tested(rng_person, t_morning, test_params);
    EXPECT_EQ(person.get_infection_state(t_morning), mio::abm::InfectionState::InfectedSymptoms);
    EXPECT_EQ(mio::abm::go_to_work(rng_person, person, t_morning, dt, infection_parameters),
              mio::abm::LocationType::Home);
    EXPECT_EQ(person.get_infection_state(t_morning + dt), mio::abm::InfectionState::Recovered);

    // Test removal from quarantine.
    person.remove_quarantine();
    EXPECT_EQ(mio::abm::go_to_work(rng_person, person, t_morning, dt, infection_parameters),
              mio::abm::LocationType::Work);
}

/**
 * @brief Test the get_tested function for both infected and susceptible individuals.
 */
TEST_F(TestPerson, get_tested)
{
    using testing::Return;
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);

    mio::abm::TimePoint t(0);
    mio::abm::Location loc(mio::abm::LocationType::Home, 0, num_age_groups);
    auto infected =
        make_test_person(this->get_rng(), loc, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms);
    auto rng_infected   = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), infected);
    auto susceptible    = mio::abm::Person(this->get_rng(), loc.get_type(), loc.get_id(), age_group_15_to_34);
    auto rng_suscetible = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), susceptible);

    auto pcr_parameters     = params.get<mio::abm::TestData>()[mio::abm::TestType::PCR];
    auto antigen_parameters = params.get<mio::abm::TestData>()[mio::abm::TestType::Antigen];
    // Test pcr test
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>>
        mock_uniform_dist_pcr;
    EXPECT_CALL(mock_uniform_dist_pcr.get_mock(), invoke)
        .Times(6)
        .WillOnce(Return(0.4)) // Draw for agent's test true positive
        .WillOnce(Return(0.8)) // Draw for is_compliant() return true
        .WillOnce(Return(0.95)) // Draw for agent's test false negative
        .WillOnce(Return(0.6)) // Draw for agent's test true negative
        .WillOnce(Return(0.999)) // Draw for agent's test false negative
        .WillOnce(Return(0.8)); // Draw for is_compliant() return true

    // Verify that the infected person tests positive and is quarantined.
    EXPECT_EQ(infected.get_tested(rng_infected, t, pcr_parameters), true);
    EXPECT_EQ(infected.is_in_quarantine(t, params), true);

    // Verify that the infected person get test false negative and is not quarantined.
    infected.remove_quarantine();
    EXPECT_EQ(infected.get_tested(rng_infected, t, pcr_parameters), false);
    EXPECT_EQ(infected.is_in_quarantine(t, params), false);

    // Verify that the susceptible person tests true negative and is not quarantined.
    EXPECT_EQ(susceptible.get_tested(rng_suscetible, t, pcr_parameters), false);
    EXPECT_EQ(susceptible.is_in_quarantine(t, params), false);

    // Verify that the susceptible person tests false negative and is quarantined.
    EXPECT_EQ(susceptible.get_tested(rng_suscetible, t, pcr_parameters), true);
    EXPECT_EQ(susceptible.is_in_quarantine(t, params), true);

    // Test antigen test
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>>
        mock_uniform_dist_antigen;
    EXPECT_CALL(mock_uniform_dist_antigen.get_mock(), invoke)
        .Times(6)
        .WillOnce(Return(0.4)) // Draw for agent's test true positive
        .WillOnce(Return(0.8)) // Draw for is_compliant() return true
        .WillOnce(Return(0.95)) // Draw for agent's test false negative
        .WillOnce(Return(0.6)) // Draw for agent's test true negative
        .WillOnce(Return(0.999)) // Draw for agent's test false negative
        .WillOnce(Return(0.8)); // Draw for is_compliant() return true

    // Verify that the infected and susceptible persons get the according results.
    EXPECT_EQ(infected.get_tested(rng_infected, t, antigen_parameters), true);
    EXPECT_EQ(infected.get_tested(rng_infected, t, antigen_parameters), false);
    EXPECT_EQ(susceptible.get_tested(rng_suscetible, t, antigen_parameters), false);
    EXPECT_EQ(susceptible.get_tested(rng_suscetible, t, antigen_parameters), true);
}

/**
 * @brief Test that a Person can change locations and correctly update cell indices.
 */
TEST_F(TestPerson, getCells)
{
    // Initialize home and target locations with designated cells.
    mio::abm::Location home(mio::abm::LocationType::Home, 0, 6, 1);
    mio::abm::Location location(mio::abm::LocationType::PublicTransport, 1, 6, 7);
    // Create a test person at the home location.
    auto person =
        make_test_person(this->get_rng(), home, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);

    // Move the person to a new location with specified cells (3, 5).
    EXPECT_TRUE(mio::abm::change_location(person, location, mio::abm::TransportMode::Unknown, {3, 5}));

    // Check that the person's cell indices have been updated correctly.
    EXPECT_EQ(person.get_cells().size(), 2);
    EXPECT_EQ(person.get_cells()[0], 3u);
    EXPECT_EQ(person.get_cells()[1], 5u);
}

/**
 * @brief Test the interaction of a Person at a location to ensure they accumulate time correctly.
 */
TEST_F(TestPerson, interact)
{
    // Location.interact is tested seperately in the location
    auto infection_parameters = mio::abm::Parameters(num_age_groups);
    // Create a location and parameters for interaction testing.
    mio::abm::Location loc(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::TimePoint t(0);
    // Create a person and set up a random number generator specific to that person.
    auto person     = mio::abm::Person(this->get_rng(), loc.get_type(), loc.get_id(), age_group_15_to_34);
    auto rng_person = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), person);
    auto dt         = mio::abm::seconds(8640); //0.1 days
    // Simulate interaction and check that the person accumulates time at the location.
    interact_testing(rng_person, person, loc, {person}, t, dt, infection_parameters);
    EXPECT_EQ(person.get_time_at_location(), dt);
}

/**
 * @brief Test that a Person can set and verify their mask type.
 */
TEST_F(TestPerson, setWearMask)
{
    auto t = mio::abm::TimePoint(0);
    mio::abm::Location location(mio::abm::LocationType::School, 0, num_age_groups);
    auto person = make_test_person(this->get_rng(), location);

    // Test setting and verifying different mask types.
    person.set_mask(mio::abm::MaskType::None, t);
    EXPECT_EQ(person.get_mask().get_type(), mio::abm::MaskType::None);

    person.set_mask(mio::abm::MaskType::Community, t);
    EXPECT_NE(person.get_mask().get_type(), mio::abm::MaskType::None);
}

/**
 * @brief Test the mask protective factor retrieval based on model parameters.
 */
TEST_F(TestPerson, getMaskProtectiveFactor)
{
    auto t = mio::abm::TimePoint(0);
    mio::abm::Location location(mio::abm::LocationType::School, 0, 6);
    auto person_community = make_test_person(this->get_rng(), location);
    person_community.set_mask(mio::abm::MaskType::Community, t);
    auto person_surgical = make_test_person(this->get_rng(), location);
    person_surgical.set_mask(mio::abm::MaskType::Surgical, t);
    auto person_ffp2 = make_test_person(this->get_rng(), location);
    person_ffp2.set_mask(mio::abm::MaskType::FFP2, t);
    auto person_without = make_test_person(this->get_rng(), location);
    person_without.set_mask(mio::abm::MaskType::None, t);

    mio::abm::Parameters params                                             = mio::abm::Parameters(num_age_groups);
    params.get<mio::abm::MaskProtection>()[{mio::abm::MaskType::Community}] = 0.5;
    params.get<mio::abm::MaskProtection>()[{mio::abm::MaskType::Surgical}]  = 0.8;
    params.get<mio::abm::MaskProtection>()[{mio::abm::MaskType::FFP2}]      = 0.9;

    // Verify that the correct mask protection factor is returned.
    EXPECT_EQ(person_community.get_mask_protective_factor(params), 0.5);
    EXPECT_EQ(person_surgical.get_mask_protective_factor(params), 0.8);
    EXPECT_EQ(person_ffp2.get_mask_protective_factor(params), 0.9);
    EXPECT_EQ(person_without.get_mask_protective_factor(params), 0.);
}

/**
 * @brief Test retrieval of the most recent protection event for a Person.
 */
TEST_F(TestPerson, getLatestProtection)
{
    auto location               = mio::abm::Location(mio::abm::LocationType::School, 0, num_age_groups);
    auto person = mio::abm::Person(this->get_rng(), location.get_type(), location.get_id(), age_group_15_to_34);
    auto prng   = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), person);
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);

    auto t = mio::abm::TimePoint(0);
    person.add_new_vaccination(mio::abm::ProtectionType::GenericVaccine, t);
    // Verify that the latest protection is a vaccination.
    auto latest_protection = person.get_latest_protection();
    EXPECT_EQ(latest_protection.type, mio::abm::ProtectionType::GenericVaccine);
    EXPECT_EQ(latest_protection.time.days(), t.days());

    t = mio::abm::TimePoint(40 * 24 * 60 * 60);
    person.add_new_infection(mio::abm::Infection(prng, static_cast<mio::abm::VirusVariant>(0), age_group_15_to_34,
                                                 params, t, mio::abm::InfectionState::Exposed));
    latest_protection = person.get_latest_protection();
    // Verify that the latest protection is a natural infection.
    EXPECT_EQ(latest_protection.type, mio::abm::ProtectionType::NaturalInfection);
    EXPECT_EQ(latest_protection.time.days(), t.days());
}

/**
 * @brief Test that a person's RNG counter increments correctly.
 */
TEST_F(TestPerson, rng)
{
    auto p =
        mio::abm::Person(this->get_rng(), mio::abm::LocationType::Home, 0, age_group_35_to_59, mio::abm::PersonId(13));

    EXPECT_EQ(p.get_rng_counter(), mio::Counter<uint32_t>(0));

    // Verify RNG counter increments.
    auto p_rng = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), p);
    EXPECT_EQ(p_rng.get_counter(), mio::rng_totalsequence_counter<uint64_t>(13, mio::Counter<uint32_t>{0}));

    p_rng();
    EXPECT_EQ(p.get_rng_counter(), mio::Counter<uint32_t>(1));
    EXPECT_EQ(p_rng.get_counter(), mio::rng_totalsequence_counter<uint64_t>(13, mio::Counter<uint32_t>{1}));
}

/**
 * @brief Test adding and retrieving test results for a Person.
 */
TEST_F(TestPerson, addAndGetTestResult)
{
    mio::abm::Location location(mio::abm::LocationType::School, 0, num_age_groups);
    auto person = make_test_person(this->get_rng(), location);
    auto t      = mio::abm::TimePoint(0);
    // Tests if m_test_results initialized correctly
    EXPECT_EQ(person.get_test_result(mio::abm::TestType::Generic).time_of_testing,
              mio::abm::TimePoint(std::numeric_limits<int>::min()));
    EXPECT_FALSE(person.get_test_result(mio::abm::TestType::Generic).result);
    EXPECT_EQ(person.get_test_result(mio::abm::TestType::Antigen).time_of_testing,
              mio::abm::TimePoint(std::numeric_limits<int>::min()));
    EXPECT_FALSE(person.get_test_result(mio::abm::TestType::Antigen).result);
    EXPECT_EQ(person.get_test_result(mio::abm::TestType::PCR).time_of_testing,
              mio::abm::TimePoint(std::numeric_limits<int>::min()));
    EXPECT_FALSE(person.get_test_result(mio::abm::TestType::PCR).result);
    // Test if m_test_results updated
    person.add_test_result(t, mio::abm::TestType::Generic, true);
    EXPECT_TRUE(person.get_test_result(mio::abm::TestType::Generic).result);
}

/**
 * @brief Test if a Person complies with an intervention based on their compliance level.
 */
TEST_F(TestPerson, isCompliant)
{
    using testing::Return;

    // Create locations
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);

    // Create test person and associated random number generator
    auto person     = make_test_person(this->get_rng(), home);
    auto rng_person = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), person);

    // Test cases with a complete truth table for compliance levels
    struct TestCase {
        mio::abm::InterventionType intervention_type;
        double compliance_level;
        bool expected_compliance;
    };

    std::vector<TestCase> test_cases = {
        {mio::abm::InterventionType::Mask, 1.0, true},       {mio::abm::InterventionType::Mask, 0.4, false},
        {mio::abm::InterventionType::Mask, 0.2, false},      {mio::abm::InterventionType::Mask, 0.9, true},
        {mio::abm::InterventionType::Mask, 0.5, false},      {mio::abm::InterventionType::Mask, 0.1, false},
        {mio::abm::InterventionType::Testing, 1.0, true},    {mio::abm::InterventionType::Testing, 0.4, false},
        {mio::abm::InterventionType::Testing, 0.2, false},   {mio::abm::InterventionType::Testing, 0.9, true},
        {mio::abm::InterventionType::Testing, 0.5, false},   {mio::abm::InterventionType::Testing, 0.1, false},
        {mio::abm::InterventionType::Isolation, 1.0, true},  {mio::abm::InterventionType::Isolation, 0.4, false},
        {mio::abm::InterventionType::Isolation, 0.2, false}, {mio::abm::InterventionType::Isolation, 0.9, true},
        {mio::abm::InterventionType::Isolation, 0.5, false}, {mio::abm::InterventionType::Isolation, 0.1, false},
    };

    // Return mock values for all tests.
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke).Times(18).WillRepeatedly(Return(0.8));

    for (const auto& test_case : test_cases) {
        // Set the compliance level for the person
        person.set_compliance(test_case.intervention_type, test_case.compliance_level);

        // Check if the person is compliant
        EXPECT_EQ(person.is_compliant(rng_person, test_case.intervention_type), test_case.expected_compliance);
    }
}
