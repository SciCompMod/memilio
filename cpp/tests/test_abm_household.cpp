/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Sascha Korf, Khoa Nguyen
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
#include "abm/household.h"
#include "abm_helpers.h"
#include <gtest/gtest.h>

TEST(TestHouseholds, test_add_household_to_model)
{
    auto member1 = mio::abm::HouseholdMember(num_age_groups);
    member1.set_age_weight(age_group_0_to_4, 1);

    auto member2 = mio::abm::HouseholdMember(num_age_groups);
    member2.set_age_weight(age_group_5_to_14, 1);

    auto household = mio::abm::Household();
    household.add_members(member1, 2);
    household.add_members(member2, 2);

    auto model = mio::abm::Model(num_age_groups);

    add_household_to_model(model, household);
    auto persons = model.get_persons();

    // Test size
    EXPECT_EQ(persons.size(), 4);

    // Test age
    EXPECT_EQ(persons[0].get_age(), age_group_0_to_4);
    EXPECT_EQ(persons[1].get_age(), age_group_0_to_4);
    EXPECT_EQ(persons[2].get_age(), age_group_5_to_14);
    EXPECT_EQ(persons[3].get_age(), age_group_5_to_14);

    // Test location
    EXPECT_EQ(persons[0].get_location(), persons[1].get_location());
    EXPECT_EQ(persons[2].get_location(), persons[3].get_location());
}

TEST(TestHouseholds, test_add_household_group_to_model)
{

    auto member1 = mio::abm::HouseholdMember(num_age_groups);
    member1.set_age_weight(age_group_35_to_59, 1);

    auto member2 = mio::abm::HouseholdMember(num_age_groups);
    member2.set_age_weight(age_group_5_to_14, 1);

    auto household_group = mio::abm::HouseholdGroup();

    auto household1 = mio::abm::Household();
    household1.add_members(member1, 10);
    household1.add_members(member2, 2);
    household_group.add_households(household1, 5);

    auto household2 = mio::abm::Household();
    household2.add_members(member1, 2);
    household2.add_members(member2, 2);
    household_group.add_households(household2, 10);

    auto model = mio::abm::Model(num_age_groups);

    add_household_group_to_model(model, household_group);
    auto persons = model.get_persons();

    // Test size
    EXPECT_EQ(persons.size(), 100);

    // Test age
    int number_of_age5to14_year_olds = 0, number_of_age35to59_year_olds = 0;

    for (auto& person : persons) {
        if (person.get_age() == age_group_5_to_14) {
            number_of_age5to14_year_olds++;
        }
        if (person.get_age() == age_group_35_to_59) {
            number_of_age35to59_year_olds++;
        }
    }
    EXPECT_EQ(number_of_age5to14_year_olds, 30);
    EXPECT_EQ(number_of_age35to59_year_olds, 70);

    // Test location for some people
    EXPECT_EQ(persons[0].get_location(), persons[1].get_location());
    EXPECT_EQ(persons[1].get_location(), persons[5].get_location());
    EXPECT_EQ(persons[5].get_location(), persons[10].get_location());

    EXPECT_EQ(persons[60].get_location(), persons[61].get_location());
    EXPECT_EQ(persons[61].get_location(), persons[62].get_location());
    EXPECT_EQ(persons[62].get_location(), persons[63].get_location());
}
