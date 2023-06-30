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

TEST(TestPerson, init)
{
    auto location = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto t        = mio::abm::TimePoint(0);
    auto person   = mio::abm::Person(location, mio::abm::AgeGroup::Age60to79);

    ASSERT_EQ(person.get_infection_state(t), mio::abm::InfectionState::Susceptible);
    ASSERT_EQ(person.get_location(), location);
    ASSERT_EQ(person.get_person_id(), mio::abm::INVALID_PERSON_ID);
}

TEST(TestPerson, migrate)
{
    auto t      = mio::abm::TimePoint(0);
    auto home   = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto loc1   = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 1);
    auto loc2   = mio::abm::Location(mio::abm::LocationType::School, 0);
    auto loc3   = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 2);
    auto person = make_test_person(home, mio::abm::AgeGroup::Age0to4, mio::abm::InfectionState::Recovered);
    person.migrate_to(loc1, {0});

    ASSERT_EQ(person.get_location(), loc1);
    ASSERT_EQ(loc1.get_subpopulation(t, mio::abm::InfectionState::Recovered), 1);
    ASSERT_EQ(home.get_subpopulation(t, mio::abm::InfectionState::Recovered), 0);
    ASSERT_EQ(loc1.get_cells()[0].m_persons.size(), 1u);

    person.migrate_to(loc2);

    ASSERT_EQ(person.get_location(), loc2);
    ASSERT_EQ(loc2.get_subpopulation(t, mio::abm::InfectionState::Recovered), 1);
    ASSERT_EQ(loc1.get_subpopulation(t, mio::abm::InfectionState::Recovered), 0);
    ASSERT_EQ(loc1.get_cells()[0].m_persons.size(), 0u);

    person.migrate_to(loc3, {0, 1});

    ASSERT_EQ(loc3.get_cells()[0].m_persons.size(), 1u);
    ASSERT_EQ(loc3.get_cells()[1].m_persons.size(), 1u);
    ASSERT_EQ(person.get_cells().size(), 2);
    ASSERT_EQ(person.get_cells()[0], 0u);
    ASSERT_EQ(person.get_cells()[1], 1u);
}

TEST(TestPerson, setGetAssignedLocation)
{
    auto location = mio::abm::Location(mio::abm::LocationType::Work, 2);
    auto person   = mio::abm::Person(location, mio::abm::AgeGroup::Age35to59);
    person.set_assigned_location(location);
    ASSERT_EQ((int)person.get_assigned_location_index(mio::abm::LocationType::Work), 2);

    person.set_assigned_location({4, mio::abm::LocationType::Work});
    ASSERT_EQ((int)person.get_assigned_location_index(mio::abm::LocationType::Work), 4);
}

TEST(TestPerson, quarantine)
{
    using testing::Return;

    auto infection_parameters = mio::abm::GlobalInfectionParameters();
    auto home                 = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto work                 = mio::abm::Location(mio::abm::LocationType::Work, 0);

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
    infection_parameters.get<mio::abm::InfectedSymptomsToRecovered>()[{
        mio::abm::VirusVariant::Wildtype, mio::abm::AgeGroup::Age35to59, mio::abm::VaccinationState::Unvaccinated}] =
        0.5 * dt.days();

    auto person = make_test_person(home, mio::abm::AgeGroup::Age35to59, mio::abm::InfectionState::InfectedSymptoms, t_morning,
                                   infection_parameters);

    person.detect_infection(t_morning);

    ASSERT_EQ(person.get_infection_state(t_morning), mio::abm::InfectionState::InfectedSymptoms);
    ASSERT_EQ(mio::abm::go_to_work(person, t_morning, dt, {}), mio::abm::LocationType::Home);
    ASSERT_EQ(person.get_infection_state(t_morning + dt), mio::abm::InfectionState::Recovered);
    person.remove_quarantine();
    ASSERT_EQ(mio::abm::go_to_work(person, t_morning, dt, {}), mio::abm::LocationType::Work);
}

TEST(TestPerson, get_tested)
{
    using testing::Return;

    mio::abm::TimePoint t(0);
    auto loc         = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto infected    = make_test_person(loc, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedSymptoms);
    auto susceptible = mio::abm::Person(loc, mio::abm::AgeGroup::Age15to34);

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
    ASSERT_EQ(infected.get_tested(t, pcr_test.get_default()), true);
    ASSERT_EQ(infected.is_in_quarantine(), true);
    ASSERT_EQ(infected.get_tested(t, pcr_test.get_default()), false);
    ASSERT_EQ(infected.is_in_quarantine(), false);
    ASSERT_EQ(susceptible.get_tested(t, pcr_test.get_default()), false);
    ASSERT_EQ(susceptible.is_in_quarantine(), false);
    ASSERT_EQ(susceptible.get_tested(t, pcr_test.get_default()), true);
    ASSERT_EQ(susceptible.is_in_quarantine(), true);
    ASSERT_EQ(susceptible.get_time_since_negative_test(), mio::abm::days(0));

    // Test antigen test
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>>
        mock_uniform_dist_antigen;
    EXPECT_CALL(mock_uniform_dist_antigen.get_mock(), invoke)
        .Times(4)
        .WillOnce(Return(0.4))
        .WillOnce(Return(0.95))
        .WillOnce(Return(0.6))
        .WillOnce(Return(0.999));
    ASSERT_EQ(infected.get_tested(t, antigen_test.get_default()), true);
    ASSERT_EQ(infected.get_tested(t, antigen_test.get_default()), false);
    ASSERT_EQ(susceptible.get_tested(t, antigen_test.get_default()), false);
    ASSERT_EQ(susceptible.get_tested(t, antigen_test.get_default()), true);
    ASSERT_EQ(susceptible.get_time_since_negative_test(), mio::abm::days(0));
}

TEST(TestPerson, getCells)
{
    auto home     = mio::abm::Location(mio::abm::LocationType::Home, 0, 1);
    auto location = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 2);
    auto person   = make_test_person(home, mio::abm::AgeGroup::Age15to34, mio::abm::InfectionState::InfectedNoSymptoms);
    home.add_person(person);
    person.migrate_to(location, {0, 1});
    ASSERT_EQ(person.get_cells().size(), 2);
}

TEST(TestPerson, interact)
{
    using testing::Return;
    // Location.interact is tested seperately in the location
    auto infection_parameters = mio::abm::GlobalInfectionParameters();
    auto loc                  = mio::abm::Location(mio::abm::LocationType::Home, 0);
    mio::abm::TimePoint t(0);
    auto person = mio::abm::Person(loc, mio::abm::AgeGroup::Age15to34);
    auto dt     = mio::abm::seconds(8640); //0.1 days
    person.interact(t, dt, infection_parameters);
    EXPECT_EQ(person.get_time_at_location(), dt);
}

TEST(TestPerson, applyMaskIntervention)
{
    auto home   = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto target = mio::abm::Location(mio::abm::LocationType::Work, 0);
    auto person = make_test_person(home);
    person.get_mask().change_mask(mio::abm::MaskType::Community);

    target.set_npi_active(false);
    person.apply_mask_intervention(target);
    ASSERT_FALSE(person.get_wear_mask());

    auto preferences = std::vector<double>((uint32_t)mio::abm::LocationType::Count, 1.);
    person.set_mask_preferences(preferences);
    person.apply_mask_intervention(target);

    ASSERT_TRUE(person.get_wear_mask());

    target.set_npi_active(true);
    target.set_required_mask(mio::abm::MaskType::Surgical);
    preferences = std::vector<double>((uint32_t)mio::abm::LocationType::Count, 0.);
    person.set_mask_preferences(preferences);
    person.apply_mask_intervention(target);

    ASSERT_EQ(person.get_mask().get_type(), mio::abm::MaskType::Surgical);
    ASSERT_TRUE(person.get_wear_mask());

    preferences = std::vector<double>((uint32_t)mio::abm::LocationType::Count, -1.);
    person.set_mask_preferences(preferences);
    person.apply_mask_intervention(target);

    ASSERT_FALSE(person.get_wear_mask());
}

TEST(TestPerson, setWearMask)
{
    auto location = mio::abm::Location(mio::abm::LocationType::School, 0);
    auto person   = make_test_person(location);

    person.set_wear_mask(false);
    ASSERT_FALSE(person.get_wear_mask());

    person.set_wear_mask(true);
    ASSERT_TRUE(person.get_wear_mask());
}

TEST(TestPerson, getProtectiveFactor)
{
    auto location         = mio::abm::Location(mio::abm::LocationType::School, 0);
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

    mio::abm::GlobalInfectionParameters params;
    params.get<mio::abm::MaskProtection>()[{mio::abm::MaskType::Community}] = 0.5;
    params.get<mio::abm::MaskProtection>()[{mio::abm::MaskType::Surgical}]  = 0.8;
    params.get<mio::abm::MaskProtection>()[{mio::abm::MaskType::FFP2}]      = 0.9;

    ASSERT_EQ(person_community.get_mask_protective_factor(params), 0.5);
    ASSERT_EQ(person_surgical.get_mask_protective_factor(params), 0.8);
    ASSERT_EQ(person_ffp2.get_mask_protective_factor(params), 0.9);
    ASSERT_EQ(person_without.get_mask_protective_factor(params), 0.);
}
