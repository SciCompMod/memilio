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



#include "epidemiology/abm/age.h"
#include "epidemiology/utils/custom_index_array.h"
#include "epidemiology/abm/abm.h"
#include <numeric>
#include <vector>

namespace epi {

/**
 * Overview:
 *  A household member has an vector with weighted age distribution from which the age can be calculated.
 *  A household holds a vector with household members.
 *  A household group holds a vector with a tuple of a household and the amount of times the household is in the group.
 *  E.g. if 10 households hold a member of "child and parent" and that household would be parentAndChildHousehold the vector would contain <parentAndChildHousehold, 10>. parentAndChildHousehold would consist of a vector with the members parent and child and parent has, for example the age distribution {0,0,10,1,0,0} with resprect to the AbmAgeGroup class. The child would have, for example, {1,1,0,0,0}, meaning that the age group 0-4 and 5-14 are the only possible ages and are both equally likely.
 */


/**
 * A household member represented by a weighted age distribution.
 */
class HouseholdMember{
public:
    
    /**
     * Constructs a newhousehold member.
     */
    HouseholdMember(): m_age_weights({epi::AbmAgeGroup::Count}, 0) {}
    
    /**
     * @brief Sets the weight of an age group.
     */
    void set_age_weight(epi::AbmAgeGroup age_group, int weight)
    {
        m_age_weights[{age_group}] = weight;
    }
    
    /**
     * @brief Returns the Array with the weight of each age group.
     * @returns An CustomIndexArray with the integer weights of the age groups.
     */
    const epi::CustomIndexArray<int, epi::AbmAgeGroup>& get_age_weights() const
    {
        return m_age_weights;
    }
    
private:
    epi::CustomIndexArray<int, epi::AbmAgeGroup> m_age_weights;
};

/**
 * A household group represented by a vector with household members.
 * The household may contain multiple members of the same type.
 */
class Household {
public:
    
    /**
     * Constructs a new household.
     */
    Household()
    : m_household_member_list()
    {
        m_number_of_members = 0;
    }
    
    /**
     * @brief Returns the number of members in the household.
     * @return Integer of number of members.
     */
    int get_total_number_of_members() const
    {
        return m_number_of_members;
    }
    
    /**
     * @brief Returns the number of households in the household group.
     * @return Integer of number of households.
     */
    const std::vector<std::tuple<epi::HouseholdMember, int>>& get_members() const
    {
        return m_household_member_list;
    }
    
    /**
     * Adds a number of the same members to a household.
     * @param household_member A household member from the HouseholdMember class.
     * @param number_of_members The amount of members to be added.
     */
    void add_members(HouseholdMember household_member, int number_of_members);
    
private:
    int m_number_of_members;
    std::vector<std::tuple<epi::HouseholdMember, int>> m_household_member_list;
};



/**
 * A group of households represented by different households.
 * The group may contain multiple households of the same type.
 */

class HouseholdGroup{
public:
    
    /**
     * Constructs a new household group.
     */
    HouseholdGroup()
    : m_household_list()
    {
        m_number_of_households = 0;
    }
    
    /**
     * @brief Returns the number of households in the household group.
     * @return Integer of number of households.
     */
    int get_total_number_of_households() const
    {
        return m_number_of_households;
    }
    
    /**
     * @brief Returns a vector of tuples. It contains the household and the amount of times that household is in the group.
     * @return std::vector of tuples.
     */
    const std::vector<std::tuple<Household, int>>& get_households() const
    {
        return m_household_list;
    }
    
    /**
     * Adds a number of households of the same kind, e.g. same members, to a household group.
     * @param household A household.
     * @param number_of_households The amount of same kind households.
     */
    void add_households(Household household, int number_of_households);

private:
    int m_number_of_households;
    std::vector<std::tuple<epi::Household, int>> m_household_list;
};

/**
 * Adds a specific household to the world class.
 * @param world The world class to which the household has to be added.
 * @param household The household to add to world.
 */
void add_household_to_world(epi::World& world, const epi::Household& household);

/**
 * Adds households from a household group to the world modell.
 * @param world The world class to which the group has to be added.
 * @param household_group The household group to add.
 */
void add_household_group_to_world(epi::World& world, const epi::HouseholdGroup& household_group);



}

#endif //EPI_ABM_HOUSEHOLD_H


