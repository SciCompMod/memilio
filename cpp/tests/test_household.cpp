/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Sascha Korf 
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
#include <stdio.h>
#include "epidemiology/abm/household.h"
#include "epidemiology/abm/abm.h"
#include <gtest/gtest.h>


TEST(TestHouseholds, test_add_household_to_world)
{
    auto member1 = epi::HouseholdMember();
    member1.set_age_weight(epi::AbmAgeGroup::Age0to4, 1);
    
    auto member2 = epi::HouseholdMember();
    member2.set_age_weight(epi::AbmAgeGroup::Age5to14, 1);
    
    auto household = epi::Household();
    household.add_members(member1, 2);
    household.add_members(member2, 2);
    
    auto world = epi::World();
    
    add_household_to_world(world, household);
    auto persons = world.get_persons();
    
    // Test size
    EXPECT_EQ(persons.size(),4);
    
    // Test age
    EXPECT_EQ(persons[0].get_age(),epi::AbmAgeGroup::Age0to4);
    EXPECT_EQ(persons[1].get_age(),epi::AbmAgeGroup::Age0to4);
    EXPECT_EQ(persons[2].get_age(),epi::AbmAgeGroup::Age5to14);
    EXPECT_EQ(persons[3].get_age(),epi::AbmAgeGroup::Age5to14);
    
    // Test location
    EXPECT_EQ(persons[0].get_location_id().index, persons[1].get_location_id().index);
    EXPECT_EQ(persons[2].get_location_id().index, persons[3].get_location_id().index);
}


TEST(TestHouseholds, test_add_household_group_to_world)
{
    
    auto member1 = epi::HouseholdMember();
    member1.set_age_weight(epi::AbmAgeGroup::Age35to59, 1);
    
    auto member2 = epi::HouseholdMember();
    member2.set_age_weight(epi::AbmAgeGroup::Age5to14, 1);

    auto household_group = epi::HouseholdGroup();
    
    auto household1 = epi::Household();
    household1.add_members(member1, 10);
    household1.add_members(member2, 2);
    household_group.add_households_to_group(household1, 5);

    auto household2 = epi::Household();
    household2.add_members(member1, 2);
    household2.add_members(member2, 2);
    household_group.add_households_to_group(household2, 10);
    
    auto world = epi::World();
    
    add_household_group_to_world(world, household_group);
    auto persons = world.get_persons();
    
    // Test size
    EXPECT_EQ(persons.size(), 100);
    
    
    // Test age
    int number_of_age5to14_year_olds = 0, number_of_age35to59_year_olds = 0;
    
    for (auto& person : persons) {
        if (person.get_age() == epi::AbmAgeGroup::Age5to14) {
            number_of_age5to14_year_olds++;
        }
        if (person.get_age() == epi::AbmAgeGroup::Age35to59) {
            number_of_age35to59_year_olds++;
        }
    }
    EXPECT_EQ(number_of_age5to14_year_olds, 30);
    EXPECT_EQ(number_of_age35to59_year_olds, 70);
    
    // Test location for some people
    EXPECT_EQ(persons[0].get_location_id().index, persons[1].get_location_id().index);
    EXPECT_EQ(persons[1].get_location_id().index, persons[5].get_location_id().index);
    EXPECT_EQ(persons[5].get_location_id().index, persons[10].get_location_id().index);
    
    EXPECT_EQ(persons[60].get_location_id().index, persons[61].get_location_id().index);
    EXPECT_EQ(persons[61].get_location_id().index, persons[62].get_location_id().index);
    EXPECT_EQ(persons[62].get_location_id().index, persons[63].get_location_id().index);
}
