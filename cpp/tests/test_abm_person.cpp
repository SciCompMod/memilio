/* 
* Copyright (C) 2020-2024 MEmilio
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
#include "abm/migration_rules.h"
#include "abm/person.h"
#include "abm/time.h"
#include "abm_helpers.h"
#include "memilio/utils/random_number_generator.h"

#include <gtest/gtest.h>

TEST(TestPerson, init)
{
    auto rng = mio::RandomNumberGenerator();

    mio::abm::Location location(mio::abm::LocationType::Work, 7, num_age_groups);
    auto t      = mio::abm::TimePoint(0);
    auto person = mio::abm::Person(rng, location.get_type(), location.get_id(), age_group_60_to_79);

    EXPECT_EQ(person.get_infection_state(t), mio::abm::InfectionState::Susceptible);
    EXPECT_EQ(person.get_location(), location.get_id());
    EXPECT_EQ(person.get_id(), mio::abm::PersonId::invalid_id());
}

TEST(TestPerson, migrate)
{
    auto rng = mio::RandomNumberGenerator();

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location loc1(mio::abm::LocationType::PublicTransport, 1, 6, 1);
    mio::abm::Location loc2(mio::abm::LocationType::School, 2, num_age_groups);
    mio::abm::Location loc3(mio::abm::LocationType::PublicTransport, 3, 6, 2);
    auto person = make_test_person(home, age_group_0_to_4, mio::abm::InfectionState::Recovered);

    // check that a person does not move to its current location
    person.add_time_at_location(mio::abm::hours(1));
    EXPECT_FALSE(mio::abm::migrate(person, home));
    EXPECT_EQ(person.get_time_at_location(), mio::abm::hours(1));
    EXPECT_EQ(person.get_location(), home.get_id());

    // move the person around a bit
    EXPECT_TRUE(mio::abm::migrate(person, loc1, mio::abm::TransportMode::Unknown, {0}));
    EXPECT_EQ(person.get_time_at_location(), mio::abm::TimeSpan(0));
    EXPECT_EQ(person.get_location(), loc1.get_id());
    EXPECT_EQ(person.get_last_transport_mode(), mio::abm::TransportMode::Unknown);

    EXPECT_TRUE(mio::abm::migrate(person, loc2, mio::abm::TransportMode::Walking, {0}));
    EXPECT_EQ(person.get_time_at_location(), mio::abm::TimeSpan(0));
    EXPECT_EQ(person.get_location(), loc2.get_id());
    EXPECT_EQ(person.get_last_transport_mode(), mio::abm::TransportMode::Walking);

    EXPECT_TRUE(mio::abm::migrate(person, loc3, mio::abm::TransportMode::Bike, {0, 1}));
    EXPECT_EQ(person.get_time_at_location(), mio::abm::TimeSpan(0));
    EXPECT_EQ(person.get_location(), loc3.get_id());
    EXPECT_EQ(person.get_last_transport_mode(), mio::abm::TransportMode::Bike);
    ASSERT_EQ(person.get_cells().size(), 2);
    EXPECT_EQ(person.get_cells()[0], 0u);
    EXPECT_EQ(person.get_cells()[1], 1u);
}

TEST(TestPerson, setGetAssignedLocation)
{
    auto rng = mio::RandomNumberGenerator();
    mio::abm::Location location(mio::abm::LocationType::Work, 2, num_age_groups);
    auto person = mio::abm::Person(rng, location.get_type(), location.get_id(), age_group_35_to_59);
    person.set_assigned_location(location.get_type(), location.get_id());
    EXPECT_EQ(person.get_assigned_location(mio::abm::LocationType::Work), mio::abm::LocationId(2));

    person.set_assigned_location(mio::abm::LocationType::Work, mio::abm::LocationId(4));
    EXPECT_EQ(person.get_assigned_location(mio::abm::LocationType::Work), mio::abm::LocationId(4));
}

TEST(TestPerson, quarantine)
{
    using testing::Return;
    auto rng         = mio::RandomNumberGenerator();
    auto test_params = mio::abm::TestParameters{1.01, 1.01}; //100% safe test

    auto infection_parameters = mio::abm::Parameters(num_age_groups);
    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location work(mio::abm::LocationType::Work, 1, num_age_groups);

    //setup rng mock so the person has a state transition to Recovered
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(4))
        .WillOnce(testing::Return(0.6)) // workgroup
        .WillOnce(testing::Return(0.6)) // schoolgroup
        .WillOnce(testing::Return(0.6)) // goto_work_hour
        .WillOnce(testing::Return(0.6)) // goto_school_hour
        .WillRepeatedly(testing::Return(1.0)); // ViralLoad draws

    auto t_morning = mio::abm::TimePoint(0) + mio::abm::hours(7);
    auto dt        = mio::abm::hours(1);
    infection_parameters
        .get<mio::abm::InfectedSymptomsToRecovered>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] =
        0.5 * dt.days();
    infection_parameters.get<mio::abm::AgeGroupGotoSchool>().set_multiple({age_group_5_to_14}, true);
    infection_parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({age_group_15_to_34, age_group_35_to_59}, true);

    auto person     = make_test_person(home, age_group_35_to_59, mio::abm::InfectionState::InfectedSymptoms, t_morning,
                                       infection_parameters);
    auto rng_person = mio::abm::PersonalRandomNumberGenerator(rng, person);

    person.get_tested(rng_person, t_morning, test_params);

    EXPECT_EQ(person.get_infection_state(t_morning), mio::abm::InfectionState::InfectedSymptoms);
    EXPECT_EQ(mio::abm::go_to_work(rng_person, person, t_morning, dt, infection_parameters),
              mio::abm::LocationType::Home);
    EXPECT_EQ(person.get_infection_state(t_morning + dt), mio::abm::InfectionState::Recovered);
    person.remove_quarantine();
    EXPECT_EQ(mio::abm::go_to_work(rng_person, person, t_morning, dt, infection_parameters),
              mio::abm::LocationType::Work);
}

TEST(TestPerson, get_tested)
{
    using testing::Return;
    auto rng                    = mio::RandomNumberGenerator();
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);

    mio::abm::TimePoint t(0);
    mio::abm::Location loc(mio::abm::LocationType::Home, 0, num_age_groups);
    auto infected       = make_test_person(loc, age_group_15_to_34, mio::abm::InfectionState::InfectedSymptoms);
    auto rng_infected   = mio::abm::PersonalRandomNumberGenerator(rng, infected);
    auto susceptible    = mio::abm::Person(rng, loc.get_type(), loc.get_id(), age_group_15_to_34);
    auto rng_suscetible = mio::abm::PersonalRandomNumberGenerator(rng, susceptible);

    auto pcr_test     = mio::abm::PCRTest();
    auto antigen_test = mio::abm::AntigenTest();

    // Test pcr test
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>>
        mock_uniform_dist_pcr;
    EXPECT_CALL(mock_uniform_dist_pcr.get_mock(), invoke)
        .Times(4)
        .WillOnce(Return(0.4))
        .WillOnce(Return(0.95))
        .WillOnce(Return(0.6))
        .WillOnce(Return(0.999));
    EXPECT_EQ(infected.get_tested(rng_infected, t, pcr_test.get_default()), true);
    EXPECT_EQ(infected.is_in_quarantine(t, params), true);
    infected.remove_quarantine();
    EXPECT_EQ(infected.get_tested(rng_infected, t, pcr_test.get_default()), false);
    EXPECT_EQ(infected.is_in_quarantine(t, params), false);
    EXPECT_EQ(susceptible.get_tested(rng_suscetible, t, pcr_test.get_default()), false);
    EXPECT_EQ(susceptible.is_in_quarantine(t, params), false);
    EXPECT_EQ(susceptible.get_tested(rng_suscetible, t, pcr_test.get_default()), true);
    EXPECT_EQ(susceptible.is_in_quarantine(t, params), true);
    EXPECT_EQ(susceptible.get_time_of_last_test(), mio::abm::TimePoint(0));

    // Test antigen test
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>>
        mock_uniform_dist_antigen;
    EXPECT_CALL(mock_uniform_dist_antigen.get_mock(), invoke)
        .Times(4)
        .WillOnce(Return(0.4))
        .WillOnce(Return(0.95))
        .WillOnce(Return(0.6))
        .WillOnce(Return(0.999));
    EXPECT_EQ(infected.get_tested(rng_infected, t, antigen_test.get_default()), true);
    EXPECT_EQ(infected.get_tested(rng_infected, t, antigen_test.get_default()), false);
    EXPECT_EQ(susceptible.get_tested(rng_suscetible, t, antigen_test.get_default()), false);
    EXPECT_EQ(susceptible.get_tested(rng_suscetible, t, antigen_test.get_default()), true);
    EXPECT_EQ(susceptible.get_time_of_last_test(), mio::abm::TimePoint(0));
}

TEST(TestPerson, getCells)
{
    mio::abm::Location home(mio::abm::LocationType::Home, 0, 6, 1);
    mio::abm::Location location(mio::abm::LocationType::PublicTransport, 1, 6, 7);
    auto person = make_test_person(home, age_group_15_to_34, mio::abm::InfectionState::InfectedNoSymptoms);

    EXPECT_TRUE(mio::abm::migrate(person, location, mio::abm::TransportMode::Unknown, {3, 5}));

    ASSERT_EQ(person.get_cells().size(), 2);
    EXPECT_EQ(person.get_cells()[0], 3u);
    EXPECT_EQ(person.get_cells()[1], 5u);
}

TEST(TestPerson, interact)
{
    auto rng = mio::RandomNumberGenerator();

    // Location.interact is tested seperately in the location
    auto infection_parameters = mio::abm::Parameters(num_age_groups);
    mio::abm::Location loc(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::TimePoint t(0);
    auto person     = mio::abm::Person(rng, loc.get_type(), loc.get_id(), age_group_15_to_34);
    auto rng_person = mio::abm::PersonalRandomNumberGenerator(rng, person);
    auto dt         = mio::abm::seconds(8640); //0.1 days
    interact_testing(rng_person, person, loc, {person}, t, dt, infection_parameters);
    EXPECT_EQ(person.get_time_at_location(), dt);
}

TEST(TestPerson, applyMaskIntervention)
{
    auto rng = mio::RandomNumberGenerator();

    mio::abm::Location home(mio::abm::LocationType::Home, 0, num_age_groups);
    mio::abm::Location target(mio::abm::LocationType::Work, 0, num_age_groups);
    auto person = make_test_person(home);
    person.get_mask().change_mask(mio::abm::MaskType::Community);
    auto rng_person = mio::abm::PersonalRandomNumberGenerator(rng, person);

    target.set_npi_active(false);
    person.apply_mask_intervention(rng_person, target);
    ASSERT_FALSE(person.get_wear_mask());

    auto preferences = std::vector<double>((uint32_t)mio::abm::LocationType::Count, 1.);
    person.set_mask_preferences(preferences);
    person.apply_mask_intervention(rng_person, target);

    ASSERT_TRUE(person.get_wear_mask());

    target.set_npi_active(true);
    target.set_required_mask(mio::abm::MaskType::Surgical);
    preferences = std::vector<double>((uint32_t)mio::abm::LocationType::Count, 0.);
    person.set_mask_preferences(preferences);
    person.apply_mask_intervention(rng_person, target);

    ASSERT_EQ(person.get_mask().get_type(), mio::abm::MaskType::Surgical);
    ASSERT_TRUE(person.get_wear_mask());

    preferences = std::vector<double>((uint32_t)mio::abm::LocationType::Count, -1.);
    person.set_mask_preferences(preferences);
    person.apply_mask_intervention(rng_person, target);

    ASSERT_FALSE(person.get_wear_mask());
}

TEST(TestPerson, setWearMask)
{
    mio::abm::Location location(mio::abm::LocationType::School, 0, num_age_groups);
    auto person = make_test_person(location);

    person.set_wear_mask(false);
    ASSERT_FALSE(person.get_wear_mask());

    person.set_wear_mask(true);
    ASSERT_TRUE(person.get_wear_mask());
}

TEST(TestPerson, getMaskProtectiveFactor)
{
    mio::abm::Location location(mio::abm::LocationType::School, 0, 6);
    auto person_community = make_test_person(location);
    person_community.get_mask().change_mask(mio::abm::MaskType::Community);
    person_community.set_wear_mask(true);
    auto person_surgical = make_test_person(location);
    person_surgical.get_mask().change_mask(mio::abm::MaskType::Surgical);
    person_surgical.set_wear_mask(true);
    auto person_ffp2 = make_test_person(location);
    person_ffp2.get_mask().change_mask(mio::abm::MaskType::FFP2);
    person_ffp2.set_wear_mask(true);
    auto person_without = make_test_person(location);
    person_without.set_wear_mask(false);

    mio::abm::Parameters params                                             = mio::abm::Parameters(num_age_groups);
    params.get<mio::abm::MaskProtection>()[{mio::abm::MaskType::Community}] = 0.5;
    params.get<mio::abm::MaskProtection>()[{mio::abm::MaskType::Surgical}]  = 0.8;
    params.get<mio::abm::MaskProtection>()[{mio::abm::MaskType::FFP2}]      = 0.9;

    ASSERT_EQ(person_community.get_mask_protective_factor(params), 0.5);
    ASSERT_EQ(person_surgical.get_mask_protective_factor(params), 0.8);
    ASSERT_EQ(person_ffp2.get_mask_protective_factor(params), 0.9);
    ASSERT_EQ(person_without.get_mask_protective_factor(params), 0.);
}

TEST(TestPerson, getLatestProtection)
{
    auto rng                    = mio::RandomNumberGenerator();
    auto location               = mio::abm::Location(mio::abm::LocationType::School, 0, num_age_groups);
    auto person                 = mio::abm::Person(rng, location.get_type(), location.get_id(), age_group_15_to_34);
    auto prng                   = mio::abm::PersonalRandomNumberGenerator(rng, person);
    mio::abm::Parameters params = mio::abm::Parameters(num_age_groups);

    auto t = mio::abm::TimePoint(0);
    person.add_new_vaccination(mio::abm::ExposureType::GenericVaccine, t);
    auto latest_protection = person.get_latest_protection();
    ASSERT_EQ(latest_protection.first, mio::abm::ExposureType::GenericVaccine);
    ASSERT_EQ(latest_protection.second.days(), t.days());

    t = mio::abm::TimePoint(40 * 24 * 60 * 60);
    person.add_new_infection(mio::abm::Infection(prng, static_cast<mio::abm::VirusVariant>(0), age_group_15_to_34,
                                                 params, t, mio::abm::InfectionState::Exposed));
    latest_protection = person.get_latest_protection();
    ASSERT_EQ(latest_protection.first, mio::abm::ExposureType::NaturalInfection);
    ASSERT_EQ(latest_protection.second.days(), t.days());
}

TEST(Person, rng)
{
    auto rng = mio::RandomNumberGenerator();
    auto p   = mio::abm::Person(rng, mio::abm::LocationType::Home, 0, age_group_35_to_59, mio::abm::PersonId(13));

    EXPECT_EQ(p.get_rng_counter(), mio::Counter<uint32_t>(0));

    auto p_rng = mio::abm::PersonalRandomNumberGenerator(rng, p);
    EXPECT_EQ(p_rng.get_counter(), mio::rng_totalsequence_counter<uint64_t>(13, mio::Counter<uint32_t>{0}));

    p_rng();
    EXPECT_EQ(p.get_rng_counter(), mio::Counter<uint32_t>(1));
    EXPECT_EQ(p_rng.get_counter(), mio::rng_totalsequence_counter<uint64_t>(13, mio::Counter<uint32_t>{1}));
}
