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
#include "abm/location_type.h"
#include "abm/mask.h"
#include "abm/mask_type.h"
#include "abm/person.h"
#include "abm/state.h"
#include "test_abm.h"
#include "gmock/gmock.h"
#include <gtest/gtest.h>

TEST(TestMasks, init)
{
    auto mask = mio::abm::Mask(mio::abm::MaskType::Count);
    ASSERT_EQ(mask.get_time_used().seconds(), 0.);
}

TEST(TestMasks, setRequiredMasks)
{
    auto home = mio::abm::Location(mio::abm::LocationType::Home, 0);
    ASSERT_EQ(home.get_required_mask(), mio::abm::MaskType::Community);

    home.set_required_mask(mio::abm::MaskType::FFP2);
    ASSERT_EQ(home.get_required_mask(), mio::abm::MaskType::FFP2);
}

TEST(TestMasks, getType)
{
    auto mask = mio::abm::Mask(mio::abm::MaskType::Community);
    auto type = mask.get_type();
    ASSERT_EQ(type, mio::abm::MaskType::Community);
}

TEST(TestMasks, increaseTimeUsed)
{
    auto mask = mio::abm::Mask(mio::abm::MaskType::Community);
    auto dt   = mio::abm::hours(2);
    mask.increase_time_used(dt);
    ASSERT_EQ(mask.get_time_used(), mio::abm::hours(2));
}

TEST(TestMasks, maskUsage)
{
    auto home   = mio::abm::Location(mio::abm::LocationType::Home, 0);
    auto target = mio::abm::Location(mio::abm::LocationType::Work, 0);
    target.set_npi_active(true);
    target.set_required_mask(mio::abm::MaskType::Surgical);

    auto person = mio::abm::Person(home, mio::abm::InfectionState::Count, mio::abm::AgeGroup::Count, {});
    person.get_mask().change_mask(mio::abm::MaskType::Community);

    person.mask_usage(target);

    ASSERT_EQ(person.get_mask().get_type(), mio::abm::MaskType::Surgical);
    ASSERT_TRUE(person.get_wear_mask());

    target.set_npi_active(false);
    person.mask_usage(target);

    ASSERT_FALSE(person.get_wear_mask());

    auto preferences = mio::CustomIndexArray<double, mio::abm::LocationType>({mio::abm::LocationType::Count}, -1.);
    person.set_mask_preferences(preferences);
    person.get_mask().change_mask(mio::abm::MaskType::Community);
    person.mask_usage(target);

    ASSERT_EQ(person.get_mask().get_type(), mio::abm::MaskType::Community);
    ASSERT_FALSE(person.get_wear_mask());

    preferences = mio::CustomIndexArray<double, mio::abm::LocationType>({mio::abm::LocationType::Count}, 1.);
    person.set_mask_preferences(preferences);
    target.set_npi_active(false);
    person.mask_usage(target);

    ASSERT_TRUE(person.get_wear_mask());
}

TEST(TestMasks, MaskProtection)
{
    mio::abm::VaccinationState vaccination_state = mio::abm::VaccinationState::Unvaccinated;
    mio::abm::GlobalInfectionParameters params;

    //setup location with some chance of exposure
    auto infection_location = mio::abm::Location(mio::abm::LocationType::School, 0);
    auto susc_person1       = mio::abm::Person(infection_location, mio::abm::InfectionState::Susceptible,
                                               mio::abm::AgeGroup::Age15to34, params, vaccination_state);
    auto susc_person2       = mio::abm::Person(infection_location, mio::abm::InfectionState::Susceptible,
                                               mio::abm::AgeGroup::Age15to34, params, vaccination_state);
    auto infected1          = mio::abm::Person(infection_location, mio::abm::InfectionState::Carrier,
                                               mio::abm::AgeGroup::Age15to34, params, vaccination_state);
    infection_location.add_person(susc_person1);
    infection_location.add_person(susc_person2);
    infection_location.add_person(infected1);

    //cache precomputed results
    auto dt = mio::abm::days(1);
    infection_location.begin_step(dt, params);
    // susc_person1 wears a mask defualt protection is 1
    susc_person1.set_wear_mask(true);
    // susc_person2 does not wear a mask
    susc_person2.set_wear_mask(false);

    //mock so every exponential distr is 1 -> everyone with a strict positiv (not zero) percantage to infect will infect, see random_events.h logic
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillOnce(testing::Return((dt.days() / 2)));

    auto susc_1_new_infection_state = infection_location.interact(susc_person1, dt, params);
    auto susc_2_new_infection_state = infection_location.interact(susc_person2, dt, params);

    // The person susc_person1 should have full protection against an infection, susc_person2 not
    ASSERT_EQ(susc_1_new_infection_state, mio::abm::InfectionState::Susceptible);
    ASSERT_EQ(susc_2_new_infection_state, mio::abm::InfectionState::Exposed);
}

TEST(TestMasks, setMasks)
{
    auto location = mio::abm::Location(mio::abm::LocationType::School, 0);
    auto person   = mio::abm::Person(location, mio::abm::InfectionState::Count, mio::abm::AgeGroup::Count, {});

    person.set_wear_mask(false);
    ASSERT_FALSE(person.get_wear_mask());

    person.set_wear_mask(true);
    ASSERT_TRUE(person.get_wear_mask());
}