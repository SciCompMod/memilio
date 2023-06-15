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
#include "memilio/io/history.h"

std::string convert_loc_id_to_string(std::tuple<mio::abm::LocationType, uint32_t> tuple_id)
{
    return std::to_string(static_cast<std::uint32_t>(std::get<0>(tuple_id))) + "_" +
           std::to_string(std::get<1>(tuple_id));
}

template <typename T>
void write_log_to_file(const T& history)
{
    auto logg = history.get_log();
    // Write the results to a file.
    auto loc_id      = std::get<1>(logg);
    auto time_points = std::get<0>(logg);
    std::string input;
    std::ofstream myfile("test_output.txt");
    myfile << "Locations as numbers:\n";
    for (auto loc_id_index = 0; loc_id_index < loc_id[0].size(); ++loc_id_index) {
        myfile << convert_loc_id_to_string(loc_id[0][loc_id_index]) << "\n";
    }
    myfile << "Timepoints:\n";

    for (int t = 0; t < time_points.size(); ++t) {
        input += std::to_string(time_points[t]) + " ";
    }
    myfile << input << "\n";

    myfile.close();
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

    struct LogTimePoint : mio::LogAlways {
        using Type = double;
        static Type log(const mio::abm::Simulation& sim)
        {
            return sim.get_time().hours();
        }
    };
    struct LogLocationIds : mio::LogOnce {
        using Type = std::vector<std::tuple<mio::abm::LocationType, uint32_t>>;
        static Type log(const mio::abm::Simulation& sim)
        {
            Type location_ids{};
            for (auto& location : sim.get_world().get_locations()) {
                location_ids.push_back(std::make_tuple(location.get_type(), location.get_index()));
            }
            return location_ids;
        }
    };

    mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds> history;

    sim.advance(tmax, history);

    write_log_to_file(history);
}
