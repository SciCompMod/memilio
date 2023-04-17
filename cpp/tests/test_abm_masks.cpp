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

TEST(TestMasks, init)
{
    auto mask = mio::abm::Mask(mio::abm::MaskType::Count);
    ASSERT_EQ(mask.get_time_used().seconds(), 0.);
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

TEST(TestMasks, changeMask)
{
    auto mask = mio::abm::Mask(mio::abm::MaskType::Community);
    mask.increase_time_used(mio::abm::hours(2));
    ASSERT_EQ(mask.get_type(), mio::abm::MaskType::Community);
    ASSERT_EQ(mask.get_time_used(), mio::abm::hours(2));

    mask.change_mask(mio::abm::MaskType::Surgical);
    ASSERT_EQ(mask.get_type(), mio::abm::MaskType::Surgical);
    ASSERT_EQ(mask.get_time_used(), mio::abm::hours(0));
}

TEST(TestMasks, maskProtection)
{
    mio::abm::VaccinationState vaccination_state = mio::abm::VaccinationState::Unvaccinated;
    mio::abm::SimulationParameters params        = mio::abm::SimulationParameters(6);

    //setup location with some chance of exposure
    auto infection_location = mio::abm::Location(mio::abm::LocationType::School, 0, 6);
    auto susc_person1       = mio::abm::Person(infection_location, mio::abm::InfectionState::Susceptible, AGE_15_TO_34,
                                               params, vaccination_state);
    auto susc_person2       = mio::abm::Person(infection_location, mio::abm::InfectionState::Susceptible, AGE_15_TO_34,
                                               params, vaccination_state);
    auto infected1 = mio::abm::Person(infection_location, mio::abm::InfectionState::Carrier, AGE_15_TO_34, params,
                                      vaccination_state);
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
