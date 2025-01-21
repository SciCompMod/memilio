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
#include "abm/person.h"
#include "abm_helpers.h"
#include "random_number_test.h"

using TestMasks = RandomNumberTest;

/**
 * @brief Test initialization of a Mask object.
 */
TEST_F(TestMasks, init)
{
    auto t    = mio::abm::TimePoint(0);
    auto mask = mio::abm::Mask(mio::abm::MaskType::Count, t);
    // Test that the mask's time used is initially zero
    EXPECT_EQ(mask.get_time_used(t).seconds(), 0.);
}

/**
 * @brief Test getting the MaskType from a Mask object.
 */
TEST_F(TestMasks, getType)
{
    auto t    = mio::abm::TimePoint(0);
    auto mask = mio::abm::Mask(mio::abm::MaskType::Community, t);
    auto type = mask.get_type();
    // Test the mask's type
    EXPECT_EQ(type, mio::abm::MaskType::Community);
}

/**
 * @brief Test changing the type of a mask and resetting the time used.
 */
TEST_F(TestMasks, changeMask)
{
    auto t    = mio::abm::TimePoint(2 * 60 * 60);
    auto mask = mio::abm::Mask(mio::abm::MaskType::Community, mio::abm::TimePoint(0));
    // Test initial mask type and time used
    EXPECT_EQ(mask.get_type(), mio::abm::MaskType::Community);
    EXPECT_EQ(mask.get_time_used(t), mio::abm::hours(2));

    // Change the mask type and reset usage time
    mask.change_mask(mio::abm::MaskType::Surgical, t);

    // Test that the mask type is updated and the time used is reset to zero
    EXPECT_EQ(mask.get_type(), mio::abm::MaskType::Surgical);
    EXPECT_EQ(mask.get_time_used(t), mio::abm::hours(0));
}

/**
 * @brief Test mask protection during person interactions.
 */
TEST_F(TestMasks, maskProtection)
{
    mio::abm::Parameters params(num_age_groups);

    // Set incubation period to two days so that newly infected person is still exposed
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::LogNormalDistribution<double>>>> mock_lognorm_dist;
    EXPECT_CALL(mock_lognorm_dist.get_mock(), invoke)
        .Times(testing::AtLeast(2))
        .WillRepeatedly(testing::Return(2)); // Time in every state is two days

    // Setup location and persons for the test
    auto t = mio::abm::TimePoint(0);
    mio::abm::Location infection_location(mio::abm::LocationType::School, 0, num_age_groups);
    // Two susceptible persons, one with a mask, one without
    auto susc_person1 = mio::abm::Person(this->get_rng(), infection_location.get_type(), infection_location.get_id(),
                                         age_group_15_to_34);
    auto susc_person2 = mio::abm::Person(this->get_rng(), infection_location.get_type(), infection_location.get_id(),
                                         age_group_15_to_34);
    // Infected person interacting with susceptible persons
    auto infected1 = make_test_person(this->get_rng(), infection_location, age_group_15_to_34,
                                      mio::abm::InfectionState::InfectedSymptoms, t, params); // infected 7 days prior

    // Cache precomputed results for 1-day intervals
    auto dt = mio::abm::days(1);
    // susc_person1 wears a mask, default protection is 1
    susc_person1.set_mask(mio::abm::MaskType::FFP2, t);
    // susc_person2 does not wear a mask
    susc_person2.set_mask(mio::abm::MaskType::None, t);

    // Mock so person 2 will get infected and person 1 will not
    ScopedMockDistribution<testing::StrictMock<MockDistribution<mio::ExponentialDistribution<double>>>>
        mock_exponential_dist;
    // Person 1 interaction with full protection
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillOnce(testing::Return(1));
    auto p1_rng = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), susc_person1);
    interact_testing(p1_rng, susc_person1, infection_location, {susc_person1, susc_person2, infected1}, t, dt, params);
    // Person 2 interaction without protection
    EXPECT_CALL(mock_exponential_dist.get_mock(), invoke).WillOnce(testing::Return(0.5));
    auto p2_rng = mio::abm::PersonalRandomNumberGenerator(this->get_rng(), susc_person2);
    interact_testing(p2_rng, susc_person2, infection_location, {susc_person1, susc_person2, infected1}, t, dt, params);

    // susc_person1 (with mask) should remain susceptible
    EXPECT_EQ(susc_person1.get_infection_state(t + dt), mio::abm::InfectionState::Susceptible);
    // susc_person2 (without mask) should be exposed
    EXPECT_EQ(susc_person2.get_infection_state(t + dt), mio::abm::InfectionState::Exposed);
}
