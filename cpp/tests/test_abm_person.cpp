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

TEST(TestPerson, init)
{
    using testing::Return;

    auto location = mio::abm::Location(mio::abm::LocationType::Work, 0, 6);
    //setup rng mock so the time_until_carrier is 1.
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformIntDistribution<int>>>>
        mock_uniform_int_dist;
    EXPECT_CALL(mock_uniform_int_dist.get_mock(), invoke).Times(testing::AtLeast(1)).WillRepeatedly(Return(0));

    auto person = mio::abm::Person(location, mio::abm::InfectionState::Exposed, mio::AgeGroup(4),
                                   mio::abm::SimulationParameters(6));
    location.add_person(person);

    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Exposed);
    ASSERT_EQ(person.get_location_id().index, location.get_index());
    ASSERT_EQ(person.get_location_id().type, location.get_type());
    ASSERT_EQ(person.get_person_id(), mio::abm::INVALID_PERSON_ID);

    auto person2 = mio::abm::Person(location, mio::abm::InfectionState::Exposed, mio::AgeGroup(4),
                                    mio::abm::SimulationParameters(6), mio::abm::VaccinationState::Unvaccinated, 0);
    ASSERT_EQ(person2.get_person_id(), 0u);

    mio::abm::TimeSpan dt = mio::abm::hours(1);
    person.interact(dt, mio::abm::SimulationParameters(6), location);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Carrier);
}

TEST(TestPerson, migrate)
{
    auto home   = mio::abm::Location(mio::abm::LocationType::Home, 0, 6, 0);
    auto loc1   = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 6, 1);
    auto loc2   = mio::abm::Location(mio::abm::LocationType::School, 0, 6, 0);
    auto loc3   = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 6, 2);
    auto person = mio::abm::Person(home, mio::abm::InfectionState::Recovered_Carrier, mio::AgeGroup(0),
                                   mio::abm::SimulationParameters(6));
    home.add_person(person);
    person.migrate_to(home, loc1, {0});

    ASSERT_EQ(person.get_location_id().index, loc1.get_index());
    ASSERT_EQ(person.get_location_id().type, loc1.get_type());
    ASSERT_EQ(loc1.get_subpopulation(mio::abm::InfectionState::Recovered_Carrier), 1);
    ASSERT_EQ(home.get_subpopulation(mio::abm::InfectionState::Recovered_Carrier), 0);
    ASSERT_EQ(loc1.get_cells()[0].num_people, 1u);

    person.migrate_to(loc1, loc2);

    ASSERT_EQ(person.get_location_id().index, loc2.get_index());
    ASSERT_EQ(person.get_location_id().type, loc2.get_type());
    ASSERT_EQ(loc2.get_subpopulation(mio::abm::InfectionState::Recovered_Carrier), 1);
    ASSERT_EQ(loc1.get_subpopulation(mio::abm::InfectionState::Recovered_Carrier), 0);
    ASSERT_EQ(loc1.get_cells()[0].num_people, 0u);

    person.migrate_to(loc2, loc3, {0, 1});

    ASSERT_EQ(loc3.get_cells()[0].num_people, 1u);
    ASSERT_EQ(loc3.get_cells()[1].num_people, 1u);
    ASSERT_EQ(person.get_cells().size(), 2);
    ASSERT_EQ(person.get_cells()[0], 0u);
    ASSERT_EQ(person.get_cells()[1], 1u);
}

TEST(TestPerson, setGetAssignedLocation)
{
    auto location = mio::abm::Location(mio::abm::LocationType::Work, 2, 6);
    auto person   = mio::abm::Person(location, mio::abm::InfectionState::Recovered_Carrier, mio::AgeGroup(4),
                                     mio::abm::SimulationParameters(6));
    person.set_assigned_location(location);
    ASSERT_EQ((int)person.get_assigned_location_index(mio::abm::LocationType::Work), 2);

    person.set_assigned_location({4, mio::abm::LocationType::Work});
    ASSERT_EQ((int)person.get_assigned_location_index(mio::abm::LocationType::Work), 4);
}

TEST(TestPerson, quarantine)
{
    using testing::Return;

    auto infection_parameters = mio::abm::SimulationParameters(6);
    auto home                 = mio::abm::Location(mio::abm::LocationType::Home, 0, 6);
    auto work                 = mio::abm::Location(mio::abm::LocationType::Work, 0, 6);

    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::UniformDistribution<double>>>> mock_uniform_dist;
    EXPECT_CALL(mock_uniform_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillOnce(testing::Return(0.6))
        .WillRepeatedly(testing::Return(1.0));

    auto person = mio::abm::Person(home, mio::abm::InfectionProperties(mio::abm::InfectionState::Infected, true),
                                   mio::AgeGroup(2), infection_parameters);
    home.add_person(person);
    auto t_morning                        = mio::abm::TimePoint(0) + mio::abm::hours(7);
    auto dt                               = mio::abm::hours(1);
    mio::abm::SimulationParameters params = mio::abm::SimulationParameters(6);
    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    params.get<mio::abm::AgeGroupGotoSchool>() = {mio::AgeGroup(1)};
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 or 35-59)
    params.get<mio::abm::AgeGroupGotoWork>() = {mio::AgeGroup(2), mio::AgeGroup(3)};
    ASSERT_EQ(mio::abm::go_to_work(person, t_morning, dt, params), mio::abm::LocationType::Home);
    //setup rng mock so the person has a state transition to Recovered_Infected
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.04));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));
    person.interact(dt, infection_parameters, home);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Recovered_Infected);
    ASSERT_EQ(mio::abm::go_to_work(person, t_morning, dt, params), mio::abm::LocationType::Work);
}

TEST(TestPerson, get_tested)
{
    using testing::Return;

    auto loc = mio::abm::Location(mio::abm::LocationType::Home, 0, 6);
    auto infected =
        mio::abm::Person(loc, mio::abm::InfectionState::Infected, mio::AgeGroup(2), mio::abm::SimulationParameters(6));
    auto susceptible = mio::abm::Person(loc, mio::abm::InfectionState::Susceptible, mio::AgeGroup(2),
                                        mio::abm::SimulationParameters(6));

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
    ASSERT_EQ(infected.get_tested(pcr_test.get_default()), true);
    ASSERT_EQ(infected.is_in_quarantine(), true);
    ASSERT_EQ(infected.get_tested(pcr_test.get_default()), false);
    ASSERT_EQ(infected.is_in_quarantine(), false);
    ASSERT_EQ(susceptible.get_tested(pcr_test.get_default()), false);
    ASSERT_EQ(susceptible.is_in_quarantine(), false);
    ASSERT_EQ(susceptible.get_tested(pcr_test.get_default()), true);
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
    ASSERT_EQ(infected.get_tested(antigen_test.get_default()), true);
    ASSERT_EQ(infected.get_tested(antigen_test.get_default()), false);
    ASSERT_EQ(susceptible.get_tested(antigen_test.get_default()), false);
    ASSERT_EQ(susceptible.get_tested(antigen_test.get_default()), true);
    ASSERT_EQ(susceptible.get_time_since_negative_test(), mio::abm::days(0));
}

TEST(TestPerson, getCells)
{
    auto home     = mio::abm::Location(mio::abm::LocationType::Home, 0, 0);
    auto location = mio::abm::Location(mio::abm::LocationType::PublicTransport, 0, 2);
    auto person =
        mio::abm::Person(home, mio::abm::InfectionState::Carrier, mio::AgeGroup(2), mio::abm::SimulationParameters(6));
    home.add_person(person);
    person.migrate_to(home, location, {0, 1});
    ASSERT_EQ(person.get_cells().size(), 2);
}

TEST(TestPerson, interact)
{
    using testing::Return;

    auto infection_parameters = mio::abm::SimulationParameters(6);
    auto loc                  = mio::abm::Location(mio::abm::LocationType::Home, 0, 6);
    auto person = mio::abm::Person(loc, mio::abm::InfectionState::Infected, mio::AgeGroup(2), infection_parameters);
    loc.add_person(person);
    auto dt = mio::abm::seconds(8640); //0.1 days
    loc.begin_step(dt, mio::abm::SimulationParameters(6));

    //setup rng mock so the person has a state transition
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.09));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    person.interact(dt, infection_parameters, loc);
    EXPECT_EQ(person.get_infection_state(), mio::abm::InfectionState::Recovered_Infected);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Recovered_Infected), 1);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Infected), 0);
}

TEST(TestPerson, interact_exposed)
{
    using testing::Return;

    auto infection_parameters = mio::abm::SimulationParameters(6);
    infection_parameters.set<mio::abm::IncubationPeriod>({{mio::AgeGroup(6), mio::abm::VaccinationState::Count}, 2.});

    //setup location with some chance of exposure
    auto loc       = mio::abm::Location(mio::abm::LocationType::Work, 0, 6);
    auto infected1 = mio::abm::Person(loc, mio::abm::InfectionState::Carrier, mio::AgeGroup(2), infection_parameters);
    loc.add_person(infected1);
    auto infected2 = mio::abm::Person(loc, mio::abm::InfectionState::Infected, mio::AgeGroup(1), infection_parameters);
    loc.add_person(infected2);
    auto infected3 = mio::abm::Person(loc, mio::abm::InfectionState::Infected, mio::AgeGroup(4), infection_parameters);
    loc.add_person(infected3);
    auto person = mio::abm::Person(loc, mio::abm::InfectionState::Susceptible, mio::AgeGroup(2), infection_parameters);
    loc.add_person(person);
    loc.begin_step(mio::abm::hours(1), mio::abm::SimulationParameters(6));

    //setup rng mock so the person becomes exposed
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::DiscreteDistribution<size_t>>>> mock_discrete_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).Times(1).WillOnce(Return(0.49));
    EXPECT_CALL(mock_discrete_dist.get_mock(), invoke).Times(1).WillOnce(Return(0));

    //person becomes exposed
    person.interact(mio::abm::hours(12), infection_parameters, loc);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Exposed);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Exposed), 1);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Carrier), 1);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Infected), 2);

    //person becomes a carrier after the incubation time runs out, not random
    person.interact(mio::abm::hours(12), infection_parameters, loc);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Exposed);

    person.interact(mio::abm::hours(12), infection_parameters, loc);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Exposed);

    person.interact(mio::abm::hours(24), infection_parameters, loc);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Exposed);

    person.interact(mio::abm::hours(1), infection_parameters, loc);
    ASSERT_EQ(person.get_infection_state(), mio::abm::InfectionState::Carrier);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Exposed), 0);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Carrier), 2);
    EXPECT_EQ(loc.get_subpopulation(mio::abm::InfectionState::Infected), 2);
}

TEST(TestPerson, applyMaskIntervention)
{
    auto home   = mio::abm::Location(mio::abm::LocationType::Home, 0, 6);
    auto target = mio::abm::Location(mio::abm::LocationType::Work, 0, 6);
    auto person =
        mio::abm::Person(home, mio::abm::InfectionState::Count, mio::AgeGroup(6), mio::abm::SimulationParameters(6));
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
    auto location = mio::abm::Location(mio::abm::LocationType::School, 0, 6);
    auto person   = mio::abm::Person(location, mio::abm::InfectionState::Count, mio::AgeGroup(6),
                                     mio::abm::SimulationParameters(6));

    person.set_wear_mask(false);
    ASSERT_FALSE(person.get_wear_mask());

    person.set_wear_mask(true);
    ASSERT_TRUE(person.get_wear_mask());
}

TEST(TestPerson, getProtectiveFactor)
{
    auto location         = mio::abm::Location(mio::abm::LocationType::School, 0, 6);
    auto person_community = mio::abm::Person(location, mio::abm::InfectionState::Count, mio::AgeGroup(6),
                                             mio::abm::SimulationParameters(6));
    person_community.get_mask().change_mask(mio::abm::MaskType::Community);
    person_community.set_wear_mask(true);
    auto person_surgical = mio::abm::Person(location, mio::abm::InfectionState::Count, mio::AgeGroup(6),
                                            mio::abm::SimulationParameters(6));
    person_surgical.get_mask().change_mask(mio::abm::MaskType::Surgical);
    person_surgical.set_wear_mask(true);
    auto person_ffp2 = mio::abm::Person(location, mio::abm::InfectionState::Count, mio::AgeGroup(6),
                                        mio::abm::SimulationParameters(6));
    person_ffp2.get_mask().change_mask(mio::abm::MaskType::FFP2);
    person_ffp2.set_wear_mask(true);
    auto person_without = mio::abm::Person(location, mio::abm::InfectionState::Count, mio::AgeGroup(6),
                                           mio::abm::SimulationParameters(6));
    person_without.set_wear_mask(false);

    mio::abm::SimulationParameters params                                   = mio::abm::SimulationParameters(6);
    params.get<mio::abm::MaskProtection>()[{mio::abm::MaskType::Community}] = 0.5;
    params.get<mio::abm::MaskProtection>()[{mio::abm::MaskType::Surgical}]  = 0.8;
    params.get<mio::abm::MaskProtection>()[{mio::abm::MaskType::FFP2}]      = 0.9;

    ASSERT_EQ(person_community.get_protective_factor(params), 0.5);
    ASSERT_EQ(person_surgical.get_protective_factor(params), 0.8);
    ASSERT_EQ(person_ffp2.get_protective_factor(params), 0.9);
    ASSERT_EQ(person_without.get_protective_factor(params), 0.);
}
