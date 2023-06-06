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
#include <cstdio>
#include "abm/world.h"
#include "memilio/io/io.h"
#include "abm/location_type.h"
#include <fstream>
#include <string>
#include <iostream>

void write_results_to_file(const mio::abm::Simulation& sim)
{
    // The results are saved in a table with 9 rows.
    // The first row is t = time, the others correspond to the number of people with a certain infection state at this time:
    // S = Susceptible, E = Exposed, I_NS = InfectedNoSymptoms, I_Sy = InfectedSymptoms, I_Sev = InfectedSevere,
    // I_Crit = InfectedCritical, R = Recovered, D = Dead
    std::ofstream myfile("abm_minimal.txt");
    myfile << "# t S E I_NS I_Sy I_Sev I_Crit R D\n";
    for (auto i = 0; i < sim.get_result().get_num_time_points(); ++i) {
        myfile << sim.get_result().get_time(i);
        auto v = sim.get_result().get_value(i);
        for (auto j = 0; j < v.size(); ++j) {
            myfile << v[j];
            if (j < v.size() - 1) {
                myfile << " ";
            }
        }
        if (i < sim.get_result().get_num_time_points() - 1) {
            myfile << "\n";
        }
    }
    myfile.close();
    std::cout << "Results written to abm_minimal.txt" << std::endl;
}

int main()
{
    // Set global infection parameters (similar to infection parameters in SECIR model) and initialize the world
    mio::abm::GlobalInfectionParameters infection_params;

    // Set same infection parameter for all age groups. For example, the incubation period is 4 days.
    infection_params.get<mio::abm::IncubationPeriod>() = 4.;

    // Create the world with infection parameters.
    auto world = mio::abm::World(infection_params);

    // There are 3 households for each household group.
    int n_households = 3;

    // For more than 1 family households we need families. These are parents and children and randoms (which are distributed like the data we have for these households).
    auto child = mio::abm::HouseholdMember(); // A child is 50/50% 0-4 or 5-14.
    child.set_age_weight(mio::abm::AgeGroup::Age0to4, 1);
    child.set_age_weight(mio::abm::AgeGroup::Age5to14, 1);

    auto parent = mio::abm::HouseholdMember(); // A parent is 50/50% 15-34 or 35-59.
    parent.set_age_weight(mio::abm::AgeGroup::Age15to34, 1);
    parent.set_age_weight(mio::abm::AgeGroup::Age35to59, 1);

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
    world.get_individualized_location(work).get_infection_parameters().set<mio::abm::MaximumContacts>(10);

    // People can get tested at work (and do this with 0.5 probability) from time point 0 to day 30.
    auto testing_min_time = mio::abm::days(1);
    auto probability      = 0.5;
    auto start_date       = mio::abm::TimePoint(0);
    auto end_date         = mio::abm::TimePoint(0) + mio::abm::days(30);
    auto test_type        = mio::abm::AntigenTest();
    auto test_at_work     = std::vector<mio::abm::LocationType>{mio::abm::LocationType::Work};
    auto testing_criteria_work =
        std::vector<mio::abm::TestingCriteria>{mio::abm::TestingCriteria({}, test_at_work, {})};
    auto testing_scheme_work =
        mio::abm::TestingScheme(testing_criteria_work, testing_min_time, start_date, end_date, test_type, probability);
    world.get_testing_strategy().add_testing_scheme(testing_scheme_work);

    // Assign infection state to each person.
    // The infection states are chosen randomly.
    auto persons = world.get_persons();
    for (auto& person : persons) {
        mio::abm::InfectionState infection_state =
            (mio::abm::InfectionState)(rand() % ((uint32_t)mio::abm::InfectionState::Count - 1));
        if (infection_state != mio::abm::InfectionState::Susceptible)
            person.add_new_infection(mio::abm::Infection(mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         world.get_global_infection_parameters(), start_date,
                                                         infection_state));
    }

    // Assign locations to the people
    for (auto& person : persons) {
        //assign shop and event
        person.set_assigned_location(event);
        person.set_assigned_location(shop);
        //assign hospital and ICU
        person.set_assigned_location(hospital);
        person.set_assigned_location(icu);
        //assign work/school to people depending on their age
        if (person.get_age() == mio::abm::AgeGroup::Age5to14) {
            person.set_assigned_location(school);
        }
        if (person.get_age() == mio::abm::AgeGroup::Age15to34 || person.get_age() == mio::abm::AgeGroup::Age35to59) {
            person.set_assigned_location(work);
        }
    }

    // During the lockdown, social events are closed for 90% of people.
    auto t_lockdown = mio::abm::TimePoint(0) + mio::abm::days(10);
    mio::abm::close_social_events(t_lockdown, 0.9, world.get_migration_parameters());

    auto t0   = mio::abm::TimePoint(0);
    auto tmax = mio::abm::TimePoint(0) + mio::abm::days(30);
    auto sim  = mio::abm::Simulation(t0, std::move(world));

    sim.advance(tmax);

    write_results_to_file(sim);
}
