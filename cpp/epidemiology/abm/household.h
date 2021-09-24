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

#ifndef household_h
#define household_h


#endif /* household_h */

#include <numeric>
#include <vector>
#include "epidemiology/abm/age.h"
#include "epidemiology/utils/custom_index_array.h"
#include "epidemiology/abm/abm.h"

namespace epi {

/**
 * Overview:
 *  A hosuehold member has an vector with weighted age distribution from which the age can be calculated.
 *  A household holds a vector with household members.
 *  A household group holds a vector with a tuple of a household and the amount of times the household is in the group.
 *  E.g. if 10 households hold a member of "child and parent" and that household would be parentAndChildHousehold the vector would contain <parentAndChildHousehold, 10>. parentAndChildHousehold would consist of a vector with the members parent and child and parent has, for example the age distribution {0,0,10,1,0,0} with resprect to the AbmAgeGroup class. The child would have, for example, {1,1,0,0,0}, meaning that the age group 0-4 and 5-14 are the only possible ages and are both equally likely.
 */


/**
 * A household member. It consists of a weighted age distribution.
 */
class HouseholdMember{
public:
    
    HouseholdMember(): m_age_dist({epi::AbmAgeGroup::Count}, 0) {}
    
    /**
     * @brief Sets the weight of an age group.
     */
    void set_age_weight(epi::AbmAgeGroup age_group, int weight)
    {
        m_age_dist[{age_group}] = weight;
    }
    
    /**
     * @brief Returns the Array with the weight of each age group.
     * @returns An CustomIndexArray with the integer weights of the age groups.
     */
    epi::CustomIndexArray<int, epi::AbmAgeGroup> get_age_dist()
    {
        return m_age_dist;
    }
    
private:
    epi::CustomIndexArray<int, epi::AbmAgeGroup> m_age_dist;
};

/**
 * A household group. It consists of a vector with household members.
 */
class Household {
public:
    
    /**
     * @brief Returns the number of members in the household.
     * @return Integer of number of members.
     */
    int get_number_of_members()
    {
        return m_number_of_members;
    }
    
    /**
     * @brief Returns the number of households in the household group.
     * @return Integer of number of households.
     */
    std::vector<HouseholdMember> get_members()
    {
        return m_household_members;
    }
    
    void add_members(HouseholdMember household_member, int number_of_members);
    
private:
    int m_number_of_members;
    std::vector<HouseholdMember> m_household_members;
};



/**
 * A household group. It consists of a vector with households.
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
    int get_number_of_households() const
    {
        return m_number_of_households;
    }
    
    /**
     * @brief Returns a vector of tuples. It contains the household and the amount of times that household is in the group.
     * @return std::vector of tuples.
     */
    std::vector<std::tuple<Household, int>> get_households() const
    {
        return m_household_list;
    }
    
    void add_households_to_group(Household household, int count);

private:
    int m_number_of_households;
    std::vector<std::tuple<epi::Household, int>> m_household_list;
};

void add_household_to_world(epi::World& world, epi::Household& household);

void add_household_group_to_world(epi::World& world, epi::HouseholdGroup& household_group);


}



