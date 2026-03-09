/*
* Copyright (C) 2020-2026 MEmilio
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
#include "abm/common_abm_loggers.h"
#include "abm/household.h"
#include "abm/lockdown_rules.h"
#include "abm/location_type.h"
#include "abm/model.h"
#include "abm/simulation.h"
#include "memilio/io/history.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/parameter_distributions.h"
#include "boost/filesystem.hpp"

#include <iostream>

int main()
{
    mio::set_log_level(mio::LogLevel::warn);

    // Define 4 age groups.
    const size_t num_age_groups   = 4;
    const auto age_group_0_to_4   = mio::AgeGroup(0);
    const auto age_group_5_to_14  = mio::AgeGroup(1);
    const auto age_group_15_to_34 = mio::AgeGroup(2);
    const auto age_group_35_to_59 = mio::AgeGroup(3);

    // Create the model.
    auto model = mio::abm::Model(num_age_groups);
    model.parameters.get<mio::abm::TimeExposedToNoSymptoms>() = mio::ParameterDistributionLogNormal(4., 1.);
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({age_group_15_to_34, age_group_35_to_59}, true);

    // Set up households.
    int n_households = 3;

    auto child = mio::abm::HouseholdMember(num_age_groups);
    child.set_age_weight(age_group_0_to_4, 1);
    child.set_age_weight(age_group_5_to_14, 1);

    auto parent = mio::abm::HouseholdMember(num_age_groups);
    parent.set_age_weight(age_group_15_to_34, 1);
    parent.set_age_weight(age_group_35_to_59, 1);

    auto twoPersonHousehold_group = mio::abm::HouseholdGroup();
    auto twoPersonHousehold_full  = mio::abm::Household();
    twoPersonHousehold_full.add_members(child, 1);
    twoPersonHousehold_full.add_members(parent, 1);
    twoPersonHousehold_group.add_households(twoPersonHousehold_full, n_households);
    add_household_group_to_model(model, twoPersonHousehold_group);

    auto threePersonHousehold_group = mio::abm::HouseholdGroup();
    auto threePersonHousehold_full  = mio::abm::Household();
    threePersonHousehold_full.add_members(child, 1);
    threePersonHousehold_full.add_members(parent, 2);
    threePersonHousehold_group.add_households(threePersonHousehold_full, n_households);
    add_household_group_to_model(model, threePersonHousehold_group);

    // Add locations.
    auto event = model.add_location(mio::abm::LocationType::SocialEvent);
    model.get_location(event).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto hospital = model.add_location(mio::abm::LocationType::Hospital);
    model.get_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto icu = model.add_location(mio::abm::LocationType::ICU);
    model.get_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto shop = model.add_location(mio::abm::LocationType::BasicsShop);
    model.get_location(shop).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto school = model.add_location(mio::abm::LocationType::School);
    model.get_location(school).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    auto work = model.add_location(mio::abm::LocationType::Work);
    model.get_location(work).get_infection_parameters().set<mio::abm::MaximumContacts>(10);

    // Set up testing at work.
    auto start_date      = mio::abm::TimePoint(0);
    auto end_date        = mio::abm::TimePoint(0) + mio::abm::days(30);
    auto test_type       = mio::abm::TestType::Antigen;
    auto test_parameters = model.parameters.get<mio::abm::TestData>()[test_type];
    auto testing_scheme  = mio::abm::TestingScheme(mio::abm::TestingCriteria(), mio::abm::days(1), start_date,
                                                   end_date, test_parameters, 0.5);
    model.get_testing_strategy().add_scheme(mio::abm::LocationType::Work, testing_scheme);

    // Randomly assign infection states.
    for (auto& person : model.get_persons()) {
        auto rng             = mio::abm::PersonalRandomNumberGenerator(person);
        auto infection_state = (mio::abm::InfectionState)(rand() % ((uint32_t)mio::abm::InfectionState::Count - 1));
        if (infection_state != mio::abm::InfectionState::Susceptible)
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model.parameters, start_date, infection_state));
    }

    // Assign locations to people.
    for (auto& person : model.get_persons()) {
        const auto pid = person.get_id();
        model.assign_location(pid, event);
        model.assign_location(pid, shop);
        model.assign_location(pid, hospital);
        model.assign_location(pid, icu);
        if (person.get_age() == age_group_5_to_14) {
            model.assign_location(pid, school);
        }
        if (person.get_age() == age_group_15_to_34 || person.get_age() == age_group_35_to_59) {
            model.assign_location(pid, work);
        }
    }

    // Run simulation.
    auto t0   = mio::abm::TimePoint(0);
    auto tmax = mio::abm::TimePoint(0) + mio::abm::days(30);
    auto sim  = mio::abm::Simulation(t0, std::move(model));

    // Set up TimeSeriesWriter-based history loggers.
    auto num_age = sim.get_model().parameters.get_num_groups();
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionStatePerAgeGroup>
        historyInfectionStatePerAgeGroup{
            Eigen::Index((size_t)mio::abm::InfectionState::Count * num_age)};
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionPerLocationTypePerAgeGroup>
        historyInfectionPerLocationTypePerAgeGroup{
            Eigen::Index((size_t)mio::abm::LocationType::Count * num_age)};

    sim.advance(tmax, historyInfectionStatePerAgeGroup, historyInfectionPerLocationTypePerAgeGroup);

    // Extract TimeSeries from histories.
    auto infection_state_per_age_group =
        std::get<0>(historyInfectionStatePerAgeGroup.get_log());
    auto infection_per_loc_type_per_age_group =
        std::get<0>(historyInfectionPerLocationTypePerAgeGroup.get_log());

    // Save results to HDF5 files.
    auto num_groups = (int)num_age;
    boost::filesystem::path result_dir("results");

    if (!mio::create_directory(result_dir.string())) {
        std::cerr << "Error creating result directory." << std::endl;
        return 1;
    }

    auto dir_infection_state = result_dir / "infection_state_per_age_group";
    if (!mio::create_directory(dir_infection_state.string())) {
        std::cerr << "Error creating infection_state_per_age_group directory." << std::endl;
        return 1;
    }
    if (!mio::save_result({infection_state_per_age_group}, {0}, num_groups,
                          (dir_infection_state / "results_run0.h5").string())) {
        std::cerr << "Error saving infection_state_per_age_group results." << std::endl;
    }

    auto dir_infection_per_loc = result_dir / "infection_per_location_type_per_age_group";
    if (!mio::create_directory(dir_infection_per_loc.string())) {
        std::cerr << "Error creating infection_per_location_type_per_age_group directory." << std::endl;
        return 1;
    }
    if (!mio::save_result({infection_per_loc_type_per_age_group}, {0}, num_groups,
                          (dir_infection_per_loc / "results_run0.h5").string())) {
        std::cerr << "Error saving infection_per_location_type_per_age_group results." << std::endl;
    }
}
