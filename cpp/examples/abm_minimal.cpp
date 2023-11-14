/*
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Khoa Nguyen
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
#include "abm/household.h"

#include <fstream>
#include <string>
#include <iostream>

int main()
{
    // This is a minimal example with children and adults < 60 year old.
    // We divided them into 4 different age groups, which are defined as follows:
    size_t NUM_AGE_GROUPS         = 4;
    const auto AGE_GROUP_0_TO_4   = mio::AgeGroup(NUM_AGE_GROUPS - 4);
    const auto AGE_GROUP_5_TO_14  = mio::AgeGroup(NUM_AGE_GROUPS - 3);
    const auto AGE_GROUP_15_TO_34 = mio::AgeGroup(NUM_AGE_GROUPS - 2);
    const auto AGE_GROUP_35_TO_59 = mio::AgeGroup(NUM_AGE_GROUPS - 1);

    // Create the world with 4 age groups.
    auto world = mio::abm::World(NUM_AGE_GROUPS);

    // Set same infection parameter for all age groups. For example, the incubation period is 4 days.
    world.parameters.get<mio::abm::IncubationPeriod>() = 4.;

    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    world.parameters.get<mio::abm::AgeGroupGotoSchool>() = {AGE_GROUP_5_TO_14};
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 and 35-59)
    world.parameters.get<mio::abm::AgeGroupGotoWork>() = {AGE_GROUP_15_TO_34, AGE_GROUP_35_TO_59};

    // Check if the parameters satisfy their contraints.
    world.parameters.check_constraints();

    // There are 3 households for each household group.
    int n_households = 10;

    // For more than 1 family households we need families. These are parents and children and randoms (which are distributed like the data we have for these households).
    auto child = mio::abm::HouseholdMember(NUM_AGE_GROUPS); // A child is 50/50% 0-4 or 5-14.
    child.set_age_weight(AGE_GROUP_0_TO_4, 1);
    child.set_age_weight(AGE_GROUP_0_TO_4, 1);

    auto parent = mio::abm::HouseholdMember(NUM_AGE_GROUPS); // A parent is 50/50% 15-34 or 35-59.
    parent.set_age_weight(AGE_GROUP_15_TO_34, 1);
    parent.set_age_weight(AGE_GROUP_35_TO_59, 1);

    // Two-person household with one parent and one child.
    auto twoPersonHousehold_group = mio::abm::HouseholdGroup();
    auto twoPersonHousehold_full  = mio::abm::Household();
    twoPersonHousehold_full.add_members(child, 1);
    twoPersonHousehold_full.add_members(parent, 1);
    twoPersonHousehold_group.add_households(twoPersonHousehold_full, n_households);
    add_household_group_to_world(world, twoPersonHousehold_group);

    // Three-person household with two parent and one child.
    auto threePersonHousehold_group = mio::abm::HouseholdGroup();
    auto threePersonHousehold_full  = mio::abm::Household();
    threePersonHousehold_full.add_members(child, 1);
    threePersonHousehold_full.add_members(parent, 2);
    threePersonHousehold_group.add_households(threePersonHousehold_full, n_households);
    add_household_group_to_world(world, threePersonHousehold_group);

    // Add one social event with 5 maximum contacts.
    // Maximum contacs limit the number of people that a person can infect while being at this location.
    auto event = world.add_location(mio::abm::LocationType::SocialEvent);
    world.get_individualized_location(event).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    // Add hospital and ICU with 5 maximum contacs.
    auto hospital = world.add_location(mio::abm::LocationType::Hospital);
    world.get_individualized_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto icu = world.add_location(mio::abm::LocationType::ICU);
    world.get_individualized_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    // Add one supermarket, maximum constacts are assumed to be 20.
    auto shop = world.add_location(mio::abm::LocationType::BasicsShop);
    world.get_individualized_location(shop).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    // At every school, the maximum contacts are 20.
    auto school = world.add_location(mio::abm::LocationType::School);
    world.get_individualized_location(school).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    // At every workplace, maximum contacts are 10.
    auto work = world.add_location(mio::abm::LocationType::Work);
    world.get_individualized_location(work).get_infection_parameters().set<mio::abm::MaximumContacts>(20);

    // Increase aerosol transmission for all locations
    world.parameters.get<mio::abm::AerosolTransmissionRates>() = 10.0;
    // Increase contact rate for all people between 15 and 34 (i.e. people meet more often in the same location)
    world.get_individualized_location(work)
        .get_infection_parameters()
        .get<mio::abm::ContactRates>()[{AGE_GROUP_15_TO_34, AGE_GROUP_15_TO_34}] = 10.0;

    // People can get tested at work (and do this with 0.5 probability) from time point 0 to day 30.
    auto testing_min_time      = mio::abm::days(1);
    auto probability           = 0.5;
    auto start_date            = mio::abm::TimePoint(0);
    auto end_date              = mio::abm::TimePoint(0) + mio::abm::days(10);
    auto test_type             = mio::abm::AntigenTest();
    auto testing_criteria_work = mio::abm::TestingCriteria();
    auto testing_scheme_work =
        mio::abm::TestingScheme(testing_criteria_work, testing_min_time, start_date, end_date, test_type, probability);
    world.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme_work);

    // Assign infection state to each person.
    // The infection states are chosen randomly with the following distribution
    std::vector<double> population_distribution{0.5, 0.3, 0.05, 0.05, 0.05, 0.05, 0.0, 0.0};
    for (auto& person : world.get_persons()) {
        mio::abm::InfectionState infection_state =
            (mio::abm::InfectionState)mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(),
                                                                                        population_distribution);
        auto rng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         world.parameters, start_date, infection_state));
        }
    }

    // Assign locations to the people
    for (auto& person : world.get_persons()) {
        //assign shop and event
        person.set_assigned_location(event);
        person.set_assigned_location(shop);
        //assign hospital and ICU
        person.set_assigned_location(hospital);
        person.set_assigned_location(icu);
        //assign work/school to people depending on their age
        if (person.get_age() == AGE_GROUP_0_TO_4) {
            person.set_assigned_location(school);
        }
        if (person.get_age() == AGE_GROUP_15_TO_34 || person.get_age() == AGE_GROUP_35_TO_59) {
            person.set_assigned_location(work);
        }
    }

    // During the lockdown, social events are closed for 90% of people.
    auto t_lockdown = mio::abm::TimePoint(0) + mio::abm::days(10);
    mio::abm::close_social_events(t_lockdown, 0.9, world.parameters);

    auto t0   = mio::abm::TimePoint(0);
    auto tmax = mio::abm::TimePoint(0) + mio::abm::days(10);
    auto sim  = mio::abm::Simulation(t0, std::move(world));

    sim.advance(tmax);

    std::ofstream outfile("abm_minimal.txt");
    sim.get_result().print_table({"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4, outfile);

    std::cout << "Results written to abm_minimal.txt" << std::endl;

    return 0;
}