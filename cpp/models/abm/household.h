/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
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

#ifndef EPI_ABM_HOUSEHOLD_H
#define EPI_ABM_HOUSEHOLD_H

#include "abm/abm.h"
#include "abm/age.h"
#include "memilio/utils/custom_index_array.h"
#include <numeric>
#include <vector>

namespace mio
{
namespace abm
{

/**
 * @file
 * A HouseholdMember has a vector with weighted age distribution from which the age can be calculated.
 * A Household holds a vector with HouseholdMembers.
 * A HouseholdGroup holds a vector with a tuple of a Household and the amount of times the Household is in the group.
 * E.g. if 10 households hold a member of "child and parent" and that Household would be parentAndChildHousehold the
 * vector would contain <parentAndChildHousehold, 10>. parentAndChildHousehold would consist of a vector with the
 * HouseholdMembers parent and child and parent has, for example the age distribution {0,0,10,1,0,0} with respect to the 
 * AgeGroup class. The child would have, for example, {1,1,0,0,0}, meaning that the AgeGroups 0-4 and 5-14 are the only
 * possible ages and are both equally likely.
 */

/**
 * @brief A HouseholdMember represented by a weighted age distribution.
 */
class HouseholdMember
{
public:
    /**
     * @brief Constructs a new HouseholdMember.
     */
    HouseholdMember()
        : m_age_weights({AgeGroup::Count}, 0)
    {
    }

    /**
     * @brief Sets the weight of an AgeGroup.
     * @todo detailed description 
     * @param[in] age_group The AgeGroup.
     * @param[in] weight The weight of the AgeGroup.
     */
    void set_age_weight(AgeGroup age_group, int weight)
    {
        m_age_weights[{age_group}] = weight;
    }

    /**
     * @brief Returns the array with the weight of each AgeGroup.
     * @returns A CustomIndexArray with the integer weights of the AgeGroups.
     */
    const CustomIndexArray<int, AgeGroup>& get_age_weights() const
    {
        return m_age_weights;
    }

private:
    CustomIndexArray<int, AgeGroup> m_age_weights;
};

/**
 * @brief A Household represented by a vector with HouseholdMembers.
 * The Household may contain multiple members of the same type.
 */
class Household
{
public:
    /**
     * @brief Constructs a new Household.
     */
    Household()
        : m_household_member_list()
    {
        m_number_of_members = 0;
        m_space_per_member  = 0;
    }

    /**
     * @brief Returns the number of members, i.e.\ Persons in the Household.
     * @return Integer of number of members.
     * @todo rework, members are not HouseholdMembers in this context right?
     */
    int get_total_number_of_members() const
    {
        return m_number_of_members;
    }

    /**
     * @brief Set the space per member that is used to compute the LocationCapacity of the Household.
     * @param[in] space_per_member space per member in cubic meters.
     * @todo rework, members are not HouseholdMembers in this context
     */
    void set_space_per_member(int space_per_member)
    {
        m_space_per_member = space_per_member;
    }

    /**
     * @brief Get the space per member of the Household.
     * @return Integer of space per member in cubic meters.
     * @todo rework, members are not HouseholdMembers in this context
     */
    int get_space_per_member() const
    {
        return m_space_per_member;
    }

    /**
     * @brief Returns the HouseholdMembers of the Household.
     * @return List of HouseholdMembers of the Household.
     */
    const std::vector<std::tuple<HouseholdMember, int>>& get_members() const
    {
        return m_household_member_list;
    }

    /**
     * @brief Adds a number of the same Householdmembers to a Household.
     * @param[in] household_member A HouseholdMember.
     * @param[in] number_of_members The amount of members to be added.
     */
    void add_members(HouseholdMember household_member, int number_of_members);

private:
    int m_number_of_members;
    int m_space_per_member; ///<space per person in cubic meters
    std::vector<std::tuple<HouseholdMember, int>> m_household_member_list;
};

/**
 * @brief A HouseholdGroup represented by different Households.
 * The group may contain multiple Households of the same type.
 */

class HouseholdGroup
{
public:
    /**
     * @brief Constructs a new HouseholdGroup.
     */
    HouseholdGroup()
        : m_household_list()
    {
        m_number_of_households = 0;
    }

    /**
     * @brief Returns the number of Households in the HouseholdGroup.
     * @return Integer of number of Households.
     */
    int get_total_number_of_households() const
    {
        return m_number_of_households;
    }

    /**
     * @brief Returns the Households of the HouseholdGroup.
     * @return a vector of tuples that contains the Houshold and the amount of times that Household is in the group.
     */
    const std::vector<std::tuple<Household, int>>& get_households() const
    {
        return m_household_list;
    }

    /**
     * @brief Adds a number of Households of the same kind, e.g.\ same members, to a HouseholdGroup.
     * @param[in] household A Household.
     * @param[in] number_of_households The amount of times that Household is in the group.
     */
    void add_households(Household household, int number_of_households);

private:
    int m_number_of_households;
    std::vector<std::tuple<Household, int>> m_household_list;
};

/**
 * @brief Adds a specific Household to the World.
 * @param world The World to which the Household has to be added.
 * @param household The Household to add to World.
 */
void add_household_to_world(World& world, const Household& household);

/**
 * @brief Adds Households from a HouseholdGroup to the World.
 * @param[out] world The World to which the group has to be added.
 * @param[in] household_group The HouseholdGroup to add.
 */
void add_household_group_to_world(World& world, const HouseholdGroup& household_group);

} // namespace abm
} // namespace mio

#endif //EPI_ABM_HOUSEHOLD_H
