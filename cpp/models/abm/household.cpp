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

#include "abm/household.h"
#include "memilio/math/eigen.h"
#include <string>

namespace mio
{
namespace abm
{

namespace
{
/**
 * @brief Picks an age from a CustomIndexArray with a weight for each AgeGroup according to a discrete distribution.
 * @param[in] age_groups A CustomIndexArray with the weights.
 * @return The picked AgeGroup.
 */
AgeGroup pick_age_group_from_age_distribution(const CustomIndexArray<int, AgeGroup>& age_groups)
{
    auto age_group_weights = age_groups.array().cast<double>().eval();
    size_t age_group       = DiscreteDistribution<size_t>::get_instance()(age_group_weights);
    return (AgeGroup)age_group;
}
} // namespace

void Household::add_members(HouseholdMember household_member, int number_of_members)
{
    m_household_member_list.push_back(std::make_tuple(household_member, number_of_members));
    m_number_of_members += number_of_members;
}

void HouseholdGroup::add_households(Household household, int number_of_households)
{
    m_household_list.push_back(std::make_tuple(household, number_of_households));
    m_number_of_households += number_of_households;
}

void add_household_to_world(World& world, const Household& household)
{
    auto home    = world.add_location(LocationType::Home);
    auto members = household.get_members();
    world.get_individualized_location(home).set_capacity(household.get_total_number_of_members(),
                                                         household.get_total_number_of_members() *
                                                             household.get_space_per_member());

    for (auto& memberTouple : members) {
        int count;
        HouseholdMember member;
        std::tie(member, count) = memberTouple;
        for (int j = 0; j < count; j++) {
            auto age_group = pick_age_group_from_age_distribution(member.get_age_weights());
            auto& person   = world.add_person(home, age_group);
            person.set_assigned_location(home);
        }
    }
}

void add_household_group_to_world(World& world, const HouseholdGroup& household_group)
{
    auto households = household_group.get_households();

    for (auto& householdTuple : households) {
        int count;
        Household household;
        std::tie(household, count) = householdTuple;
        for (int j = 0; j < count; j++) {
            add_household_to_world(world, household);
        }
    }
}

} // namespace abm
} // namespace mio
