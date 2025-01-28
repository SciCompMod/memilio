/*
* Copyright (C) 2020-2025 MEmilio
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
#include "abm/household.h"
#include "abm/lockdown_rules.h"
#include "abm/simulation.h"
#include "abm/model.h"
#include "abm/location_type.h"
#include "memilio/io/history.h"

#include <fstream>
#include <string>

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
    for (auto&& id : loc_id[0]) {
        myfile << convert_loc_id_to_string(id) << "\n";
    }
    myfile << "Timepoints:\n";

    for (auto&& t : time_points) {
        input += std::to_string(t) + " ";
    }
    myfile << input << "\n";

    myfile.close();
}

int main()
{
    // This is a minimal example with children and adults < 60y.
    // We divided them into 4 different age groups, which are defined as follows:
    const size_t num_age_groups   = 4;
    const auto age_group_0_to_4   = mio::AgeGroup(0);
    const auto age_group_5_to_14  = mio::AgeGroup(1);
    const auto age_group_15_to_34 = mio::AgeGroup(2);
    const auto age_group_35_to_59 = mio::AgeGroup(3);

    // Create the model with 4 age groups.
    auto model = mio::abm::Model(num_age_groups);

    // Set same infection parameter for all age groups. For example, the incubation period is 4 days.
    model.parameters.get<mio::abm::IncubationPeriod>() = 4.;

    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 and 35-59)
    model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({age_group_15_to_34, age_group_35_to_59}, true);

    // There are 3 households for each household group.
    int n_households = 3;

    // For more than 1 family households we need families. These are parents and children and randoms (which are distributed like the data we have for these households).
    auto child = mio::abm::HouseholdMember(num_age_groups); // A child is 50/50% 0-4 or 5-14.
    child.set_age_weight(age_group_0_to_4, 1);
    child.set_age_weight(age_group_5_to_14, 1);

    auto parent = mio::abm::HouseholdMember(num_age_groups); // A parent is 50/50% 15-34 or 35-59.
    parent.set_age_weight(age_group_15_to_34, 1);
    parent.set_age_weight(age_group_35_to_59, 1);

    // Two-person household with one parent and one child.
    auto twoPersonHousehold_group = mio::abm::HouseholdGroup();
    auto twoPersonHousehold_full  = mio::abm::Household();
    twoPersonHousehold_full.add_members(child, 1);
    twoPersonHousehold_full.add_members(parent, 1);
    twoPersonHousehold_group.add_households(twoPersonHousehold_full, n_households);
    add_household_group_to_model(model, twoPersonHousehold_group);

    // Three-person household with two parent and one child.
    auto threePersonHousehold_group = mio::abm::HouseholdGroup();
    auto threePersonHousehold_full  = mio::abm::Household();
    threePersonHousehold_full.add_members(child, 1);
    threePersonHousehold_full.add_members(parent, 2);
    threePersonHousehold_group.add_households(threePersonHousehold_full, n_households);
    add_household_group_to_model(model, threePersonHousehold_group);

    // Add one social event with 5 maximum contacts.
    // Maximum contacs limit the number of people that a person can infect while being at this location.
    auto event = model.add_location(mio::abm::LocationType::SocialEvent);
    model.get_location(event).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    // Add hospital and ICU with 5 maximum contacs.
    auto hospital = model.add_location(mio::abm::LocationType::Hospital);
    model.get_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto icu = model.add_location(mio::abm::LocationType::ICU);
    model.get_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    // Add one supermarket, maximum constacts are assumed to be 20.
    auto shop = model.add_location(mio::abm::LocationType::BasicsShop);
    model.get_location(shop).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    // At every school, the maximum contacts are 20.
    auto school = model.add_location(mio::abm::LocationType::School);
    model.get_location(school).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    // At every workplace, maximum contacts are 10.
    auto work = model.add_location(mio::abm::LocationType::Work);
    model.get_location(work).get_infection_parameters().set<mio::abm::MaximumContacts>(10);

    // People can get tested at work (and do this with 0.5 probability) from time point 0 to day 30.
    auto validity_period       = mio::abm::days(1);
    auto probability           = 0.5;
    auto start_date            = mio::abm::TimePoint(0);
    auto end_date              = mio::abm::TimePoint(0) + mio::abm::days(30);
    auto test_type             = mio::abm::TestType::Antigen;
    auto test_parameters       = model.parameters.get<mio::abm::TestData>()[test_type];
    auto testing_criteria_work = mio::abm::TestingCriteria();
    auto testing_scheme_work   = mio::abm::TestingScheme(testing_criteria_work, validity_period, start_date, end_date,
                                                         test_parameters, probability);
    model.get_testing_strategy().add_testing_scheme(mio::abm::LocationType::Work, testing_scheme_work);

    // Assign infection state to each person.
    // The infection states are chosen randomly.
    auto persons = model.get_persons();
    for (auto& person : persons) {
        auto rng = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        mio::abm::InfectionState infection_state =
            (mio::abm::InfectionState)(rand() % ((uint32_t)mio::abm::InfectionState::Count - 1));
        if (infection_state != mio::abm::InfectionState::Susceptible)
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model.parameters, start_date, infection_state));
    }

    // Assign locations to the people
    for (auto& person : model.get_persons()) {
        const auto pid = person.get_id();
        //assign shop and event
        model.assign_location(pid, event);
        model.assign_location(pid, shop);
        //assign hospital and ICU
        model.assign_location(pid, hospital);
        model.assign_location(pid, icu);
        //assign work/school to people depending on their age
        if (person.get_age() == age_group_5_to_14) {
            model.assign_location(pid, school);
        }
        if (person.get_age() == age_group_15_to_34 || person.get_age() == age_group_35_to_59) {
            model.assign_location(pid, work);
        }
    }

    // During the lockdown, social events are closed for 90% of people.
    auto t_lockdown = mio::abm::TimePoint(0) + mio::abm::days(10);
    mio::abm::close_social_events(t_lockdown, 0.9, model.parameters);

    auto t0   = mio::abm::TimePoint(0);
    auto tmax = mio::abm::TimePoint(0) + mio::abm::days(30);
    auto sim  = mio::abm::Simulation(t0, std::move(model));

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
            for (auto& location : sim.get_model().get_locations()) {
                location_ids.push_back(std::make_tuple(location.get_type(), location.get_id().get()));
            }
            return location_ids;
        }
    };

    mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds> history;

    sim.advance(tmax, history);

    write_log_to_file(history);
}
