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


#include "household.h"
#include "epidemiology/utils/eigen.h"
#include <string>


namespace epi {

/**
 * Function that adds a number of the same members to a household.
 * @param household_member A household member from the HouseholdMember class.
 * @param number_of_members The amount of members to be added.
 */
void Household::add_members(HouseholdMember household_member, int number_of_members)
{
    for (auto i = 0; i < number_of_members; i++) {
        m_household_members.push_back(household_member);
    }
    m_number_of_members += number_of_members;
}

/**
 * Function that adds a number of households of the same kind, e.g. same members, to a household group.
 * @param household A household.
 * @param count The amount of same kind households.
 */
void HouseholdGroup::add_households_to_group(Household household, int count){
    m_household_list.push_back(std::make_tuple(household, count));
    m_number_of_households += count;
}

/**
 * Function that takes a array of integers and returns a vector which has the same length with decimals. This vector has the percentage of the row with respect to the whole vector sum in each row.
 * E.g. {1,1,2} would give back {0.25,0.25,0.5}.
 * @param absoluteValueVec Array with the absolute Values.
 * @return A vector with decimals describing the percentage.
 */
std::vector<double> get_percentage_of_sum(const Eigen::Array<int, Eigen::Dynamic, 1>& absoluteValueVec){
    int sumOfValues = 0;
    std::vector<double> getPercentageOfSum;
    // Calculate sum.
    for (int i = 0; i < absoluteValueVec.size(); i++) {
        sumOfValues += absoluteValueVec(i);
    }
    // Get the percentage of each entry of the total sum.
    getPercentageOfSum.reserve(absoluteValueVec.size());
    for (int i = 0; i<absoluteValueVec.size(); i++) {
        getPercentageOfSum.push_back(((double)absoluteValueVec(i))/((double)sumOfValues));
    }
    return getPercentageOfSum;
}


/**
 * Picks an age from a custom index array with a weight for each age group according to a discrete distribution.
 * @param age_groups A custom index array with the weights.
 * @return The picked age group.
 */
epi::AbmAgeGroup pick_age_group_from_age_distribution(const epi::CustomIndexArray<int, epi::AbmAgeGroup>& age_groups){
    auto age_group_weights = age_groups.array();
    auto weights = get_percentage_of_sum(age_group_weights);
    uint32_t age_group = epi::DiscreteDistribution<size_t>::get_instance()(weights);
    return (epi::AbmAgeGroup) age_group;
}

/**
 * This function adds a specific household to the world class.
 * @param world The world class to which the household has to be added.
 * @param household The household to add to world.
 */
void add_household_to_world(epi::World& world, epi::Household& household){
    auto home = world.add_location(epi::LocationType::Home);
    auto members = household.get_members();
    for (int i = 0; i < household.get_number_of_members() ; i++) {
        auto age_group = pick_age_group_from_age_distribution(members.at(i).get_age_dist()); // Gets the age of a member from its age distribution.
        auto& person = world.add_person(home, epi::InfectionState::Susceptible, age_group); // Add person, always susceptible.
        person.set_assigned_location(home);
    }
}

/**
 * This function adds households from a household group to the world modell.
 * @param world The world class to which the group has to be added.
 * @param household_group The household group to add.
 */
void add_household_group_to_world(epi::World& world, epi::HouseholdGroup& household_group){
    auto households = household_group.get_households();
    
    for (int i = 0; i < households.size(); i++) {
        epi::Household household;
        int count;
        std::tie(household, count) = households.at(i); // Get the household (and the amount of times this household is in there) from the tuple.
        for (int j = 0; j < count; j++) {
            add_household_to_world(world, household); // Add the household count amount of times.
        }
    }
}

} // namespace epi
