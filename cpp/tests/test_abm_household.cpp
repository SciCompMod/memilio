/*
* Copyright (C) 2020-2025 MEmilio
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

/**
 * @brief Test adding a household to a model.
 * Verifies correct number of persons, age groups, and location assignments.
 */
TEST(TestHouseholds, test_add_household_to_model)
{
    // Create household members
    auto member1 = mio::abm::HouseholdMember(num_age_groups);
    member1.set_age_weight(age_group_0_to_4, 1); // Member is a child (age 0-4)

    auto member2 = mio::abm::HouseholdMember(num_age_groups);
    member2.set_age_weight(age_group_5_to_14, 1); // Member is a child (age 5-14)

    // Create household and add members
    auto household = mio::abm::Household();
    household.add_members(member1, 2); // Add two members of age group 0-4
    household.add_members(member2, 2); // Add two members of age group 5-14

    // Create model and add the household
    auto model = mio::abm::Model(num_age_groups);
    add_household_to_model(model, household);
    auto persons = model.get_persons();

    // Test: Correct number of persons added to the model
    EXPECT_EQ(persons.size(), 4);

    // Test: Correct age groups assigned to persons
    EXPECT_EQ(persons[0].get_age(), age_group_0_to_4);
    EXPECT_EQ(persons[1].get_age(), age_group_0_to_4);
    EXPECT_EQ(persons[2].get_age(), age_group_5_to_14);
    EXPECT_EQ(persons[3].get_age(), age_group_5_to_14);

    // Test: Ensure persons of the same age group are in the same location
    EXPECT_EQ(persons[0].get_location(), persons[1].get_location());
    EXPECT_EQ(persons[2].get_location(), persons[3].get_location());
}

/**
 * @brief Test adding a group of households to a model.
 * Checks correct number of persons, age distribution, and location assignments.
 */
TEST(TestHouseholds, test_add_household_group_to_model)
{
    // Create household members
    auto member1 = mio::abm::HouseholdMember(num_age_groups);
    member1.set_age_weight(age_group_35_to_59, 1); // Member is an adult (age 35-59)

    auto member2 = mio::abm::HouseholdMember(num_age_groups);
    member2.set_age_weight(age_group_5_to_14, 1); // Member is a child (age 5-14)

    // Create the first household and add members
    auto household1 = mio::abm::Household();
    household1.add_members(member1, 10); // Add ten members of age group 35-59
    household1.add_members(member2, 2); // Add two members of age group 5-14

    // Create the second household and add members
    auto household2 = mio::abm::Household();
    household2.add_members(member1, 2); // Add two members of age group 35-59
    household2.add_members(member2, 2); // Add two members of age group 5-14

    // Create a household group and add households
    auto household_group = mio::abm::HouseholdGroup();
    household_group.add_households(household1, 5); // Add household1 5 times
    household_group.add_households(household2, 10); // Add household2 10 times

    // Create model and add the household group
    auto model = mio::abm::Model(num_age_groups);
    add_household_group_to_model(model, household_group);
    auto persons = model.get_persons();

    // Test: Correct number of persons in the model (5 * 12 + 10 * 4 = 100 persons)
    EXPECT_EQ(persons.size(), 100);

    // Test: Count persons in each age group
    int number_of_age5to14_year_olds = 0, number_of_age35to59_year_olds = 0;
    for (auto& person : persons) {
        if (person.get_age() == age_group_5_to_14) {
            number_of_age5to14_year_olds++;
        }
        if (person.get_age() == age_group_35_to_59) {
            number_of_age35to59_year_olds++;
        }
    }

    // Verify the age group distribution
    EXPECT_EQ(number_of_age5to14_year_olds, 30); // 30 children (age 5-14)
    EXPECT_EQ(number_of_age35to59_year_olds, 70); // 70 adults (age 35-59)

    // Test: Ensure people in the same household share the same location (for a few checks)
    EXPECT_EQ(persons[0].get_location(), persons[1].get_location());
    EXPECT_EQ(persons[1].get_location(), persons[5].get_location());
    EXPECT_EQ(persons[5].get_location(), persons[10].get_location());

    EXPECT_EQ(persons[60].get_location(), persons[61].get_location());
    EXPECT_EQ(persons[61].get_location(), persons[62].get_location());
    EXPECT_EQ(persons[62].get_location(), persons[63].get_location());
}
