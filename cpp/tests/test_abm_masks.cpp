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
#include <gtest/gtest.h>

TEST(TestMasks, init)
{
    auto mask = mio::abm::Mask(mio::abm::MaskType::Count);
    ASSERT_EQ(mask.get_time_used().seconds(), 0.);
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
