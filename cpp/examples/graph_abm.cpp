/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Julia Bicker
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

#include "abm/abm.h"
#include <cstddef>
#include <cstdint>
#include <vector>

int main()
{
    // This is an example with three age groups representing children, adults and seniors.
    size_t num_age_groups         = 3;
    const auto age_group_children = mio::AgeGroup(0);
    const auto age_group_adults   = mio::AgeGroup(1);
    const auto age_group_seniors  = mio::AgeGroup(2);

    auto world1 = mio::abm::World(num_age_groups, 0);

    //Set infection parameters
    world1.parameters.get<mio::abm::IncubationPeriod>()              = 4.;
    world1.parameters.get<mio::abm::InfectedNoSymptomsToSymptoms>()  = 2.;
    world1.parameters.get<mio::abm::InfectedNoSymptomsToRecovered>() = 4.;
    world1.parameters.get<mio::abm::InfectedSymptomsToRecovered>()   = 5.;
    world1.parameters.get<mio::abm::InfectedSymptomsToSevere>()      = 6.;
    world1.parameters.get<mio::abm::SevereToRecovered>()             = 8.;
    world1.parameters.get<mio::abm::SevereToCritical>()              = 7.;
    world1.parameters.get<mio::abm::CriticalToRecovered>()           = 10.;
    world1.parameters.get<mio::abm::CriticalToDead>()                = 11.;

    //Age group 0 goes to school and age group 1 goes to work
    world1.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_children] = true;
    world1.parameters.get<mio::abm::AgeGroupGotoWork>()[age_group_adults]     = true;

    //Household members can be child, parent or senior
    auto child = mio::abm::HouseholdMember(num_age_groups);
    child.set_age_weight(age_group_children, 1);
    auto parent = mio::abm::HouseholdMember(num_age_groups);
    parent.set_age_weight(age_group_adults, 1);
    auto adult = mio::abm::HouseholdMember(num_age_groups);
    adult.set_age_weight(age_group_adults, 1);
    adult.set_age_weight(age_group_seniors, 1);

    //Single-Person households
    auto single_hh = mio::abm::Household();
    single_hh.add_members(adult, 1);

    //Two-Adult household
    auto two_adult_hh = mio::abm::Household();
    two_adult_hh.add_members(adult, 2);

    //Single-Parent household
    auto single_parent_hh = mio::abm::Household();
    single_parent_hh.add_members(child, 1);
    single_parent_hh.add_members(parent, 1);

    //Family household
    auto family_hh = mio::abm::Household();
    family_hh.add_members(child, 1);
    family_hh.add_members(parent, 2);

    // Vector holding all persons for the graph simulation. This vector is copied to all worlds at the end.
    std::vector<mio::abm::Person> persons;

    //Household groups for world 1
    auto single_hh_group_w1 = mio::abm::HouseholdGroup();
    single_hh_group_w1.add_households(single_hh, 5);
    auto two_adult_hh_group_w1 = mio::abm::HouseholdGroup();
    two_adult_hh_group_w1.add_households(two_adult_hh, 3);
    auto single_parent_hh_group_w1 = mio::abm::HouseholdGroup();
    single_parent_hh_group_w1.add_households(single_parent_hh, 5);
    auto family_hh_group_w1 = mio::abm::HouseholdGroup();
    family_hh_group_w1.add_households(family_hh, 10);
    add_household_group_to_world(world1, single_hh_group_w1);
    add_household_group_to_world(world1, two_adult_hh_group_w1);
    add_household_group_to_world(world1, single_hh_group_w1);
    add_household_group_to_world(world1, family_hh_group_w1);

    //add persons from world 0 to vector
    for (auto& person : world1.get_persons()) {
        persons.push_back(person);
    }

    auto world2 = mio::abm::World(num_age_groups, 1);
    //Household groups for world 2
    auto single_hh_group_w2 = mio::abm::HouseholdGroup();
    single_hh_group_w2.add_households(single_hh, 6);
    auto two_adult_hh_group_w2 = mio::abm::HouseholdGroup();
    two_adult_hh_group_w2.add_households(two_adult_hh, 2);
    auto single_parent_hh_group_w2 = mio::abm::HouseholdGroup();
    single_parent_hh_group_w2.add_households(single_parent_hh, 10);
    auto family_hh_group_w2 = mio::abm::HouseholdGroup();
    family_hh_group_w2.add_households(family_hh, 11);
    add_household_group_to_world(world2, single_hh_group_w2);
    add_household_group_to_world(world2, two_adult_hh_group_w2);
    add_household_group_to_world(world2, single_hh_group_w2);
    add_household_group_to_world(world2, family_hh_group_w2);

    //add persons from world 1 to vector
    for (auto& person : world2.get_persons()) {
        persons.push_back(person);
    }

    //Create locations for both worlds
    //world 0
    auto event_w1 = world1.add_location(mio::abm::LocationType::SocialEvent);
    world1.get_individualized_location(event_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    auto hospital_w1 = world1.add_location(mio::abm::LocationType::Hospital);
    world1.get_individualized_location(hospital_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    auto icu_w1 = world1.add_location(mio::abm::LocationType::ICU);
    world1.get_individualized_location(icu_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto shop_w1 = world1.add_location(mio::abm::LocationType::BasicsShop);
    world1.get_individualized_location(shop_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto school_w1 = world1.add_location(mio::abm::LocationType::School);
    world1.get_individualized_location(school_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto work_w1 = world1.add_location(mio::abm::LocationType::Work);
    world1.get_individualized_location(work_w1).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    //World 1
    auto event_w2 = world2.add_location(mio::abm::LocationType::SocialEvent);
    world2.get_individualized_location(event_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    auto hospital_w2 = world2.add_location(mio::abm::LocationType::Hospital);
    world2.get_individualized_location(hospital_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(10);
    auto icu_w2 = world2.add_location(mio::abm::LocationType::ICU);
    world2.get_individualized_location(icu_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto shop_w2 = world2.add_location(mio::abm::LocationType::BasicsShop);
    world2.get_individualized_location(shop_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto school_w2 = world2.add_location(mio::abm::LocationType::School);
    world2.get_individualized_location(school_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto work_w2 = world2.add_location(mio::abm::LocationType::Work);
    world2.get_individualized_location(work_w2).get_infection_parameters().set<mio::abm::MaximumContacts>(10);

    auto start_date = mio::abm::TimePoint(0);
    auto end_date   = mio::abm::TimePoint(0) + mio::abm::days(30);
    std::vector<uint32_t> params_e1;
    std::vector<uint32_t> params_e2;

    //Assign infection states and locations
    std::vector<double> infection_distribution{0.5, 0.3, 0.05, 0.05, 0.05, 0.05, 0.0, 0.0};
    for (auto& person : persons) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), infection_distribution));
        auto rng = mio::abm::Person::RandomNumberGenerator(mio::thread_local_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         world1.parameters, start_date, infection_state));
        }
        if (person.get_assigned_location_world_id(mio::abm::LocationType::Home) == world1.get_id()) {
            person.set_assigned_location(event_w1);
            person.set_assigned_location(shop_w1);
            person.set_assigned_location(hospital_w1);
            person.set_assigned_location(icu_w1);
            if (person.get_age() == age_group_children) {
                person.set_assigned_location(school_w1);
            }
            if (person.get_age() == age_group_adults) {
                //10% of adults in world 0 work in world 1
                size_t work_world = mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(),
                                                                                      std::vector<int>{90, 10});
                if (work_world == 1) { //person works in other world
                    person.set_assigned_location(work_w2);
                    //add person to edge parameters
                    params_e1.push_back(person.get_person_id());
                }
                else { //person works in same world
                    person.set_assigned_location(work_w1);
                }
            }
        }
        else {
            person.set_assigned_location(event_w2);
            person.set_assigned_location(shop_w2);
            person.set_assigned_location(hospital_w2);
            person.set_assigned_location(icu_w2);
            if (person.get_age() == age_group_children) {
                person.set_assigned_location(school_w2);
            }
            if (person.get_age() == age_group_adults) {
                //20% of adults in world 1 work in world 0
                size_t work_world = mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(),
                                                                                      std::vector<int>{20, 80});
                if (work_world == 0) { //person works in other world
                    person.set_assigned_location(work_w1);
                    //add person to edge parameters
                    params_e2.push_back(person.get_person_id());
                }
                else { //person works in same world
                    person.set_assigned_location(work_w2);
                }
            }
        }
    }

    //copy persons to both worlds
    world1.set_persons(persons);
    world2.set_persons(persons);

    return 0;
}