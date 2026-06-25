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
#include "abm/household.h"
#include "abm/lockdown_rules.h"
#include "abm/model.h"
#include "abm/common_abm_loggers.h"
#include "memilio/io/directories.h"
#include "memilio/io/history.h"

#include <fstream>

int main()
{
    // This is a minimal example with children and adults < 60 year old.
    // We divided them into 4 different age groups, which are defined as follows:
    mio::set_log_level(mio::LogLevel::warn);
    size_t num_age_groups         = 4;
    const auto age_group_0_to_4   = mio::AgeGroup(0);
    const auto age_group_5_to_14  = mio::AgeGroup(1);
    const auto age_group_15_to_34 = mio::AgeGroup(2);
    const auto age_group_35_to_59 = mio::AgeGroup(3);

    // Create the model with 4 age groups.
    auto model = mio::abm::Model(num_age_groups);
    // Set same infection parameter for all age groups. For example, the incubation period is log normally distributed with parameters 4 and 1.
    model.parameters.get<mio::abm::TimeExposedToNoSymptoms>() = mio::ParameterDistributionLogNormal(1.38, 1.);
    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 and 35-59)
    model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({age_group_15_to_34, age_group_35_to_59}, true);

    // Check if the parameters satisfy their contraints.
    model.parameters.check_constraints();

    // There are 10 households for each household group.
    int n_households = 100;

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

    // Scale location count to population size.
    // Population: ~n_households school-age kids, ~3*n_households workers, ~5*n_households total.
    // Target ~200 people per school/workplace, ~500 per shop/event.
    const int n_schools    = std::max(1, n_households / 200);
    const int n_workplaces = std::max(1, 3 * n_households / 200);
    const int n_shops      = std::max(1, 5 * n_households / 500);
    const int n_events     = std::max(1, 5 * n_households / 500);

    std::vector<mio::abm::LocationId> schools, workplaces, shops, events;
    for (int i = 0; i < n_schools; ++i) {
        auto loc = model.add_location(mio::abm::LocationType::School);
        model.get_location(loc).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
        schools.push_back(loc);
    }
    for (int i = 0; i < n_workplaces; ++i) {
        auto loc = model.add_location(mio::abm::LocationType::Work);
        model.get_location(loc).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
        model.get_location(loc).get_infection_parameters().get<mio::abm::ContactRates>().get_baseline()(
            age_group_15_to_34.get(), age_group_15_to_34.get()) = 10.0;
        workplaces.push_back(loc);
    }
    for (int i = 0; i < n_shops; ++i) {
        auto loc = model.add_location(mio::abm::LocationType::BasicsShop);
        model.get_location(loc).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
        shops.push_back(loc);
    }
    for (int i = 0; i < n_events; ++i) {
        auto loc = model.add_location(mio::abm::LocationType::SocialEvent);
        model.get_location(loc).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
        events.push_back(loc);
    }
    auto hospital = model.add_location(mio::abm::LocationType::Hospital);
    model.get_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto icu = model.add_location(mio::abm::LocationType::ICU);
    model.get_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);

    model.parameters.get<mio::abm::AerosolTransmissionRates>() = 0.10;

    // People can get tested at work (and do this with 0.5 probability) from time point 0 to day 10.
    auto validity_period       = mio::abm::days(1);
    auto probability           = 0.5;
    auto start_date            = mio::abm::TimePoint(0);
    auto end_date              = mio::abm::TimePoint(0) + mio::abm::days(10);
    auto test_type             = mio::abm::TestType::Antigen;
    auto test_parameters       = model.parameters.get<mio::abm::TestData>()[test_type];
    auto testing_criteria_work = mio::abm::TestingCriteria();
    auto testing_scheme_work   = mio::abm::TestingScheme(testing_criteria_work, validity_period, start_date, end_date,
                                                         test_parameters, probability);
    model.get_testing_strategy().add_scheme(mio::abm::LocationType::Work, testing_scheme_work);

    // Assign infection state to each person.
    std::vector<ScalarType> infection_distribution{0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

    for (auto& person : model.get_persons()) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), infection_distribution));
        auto rng = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model.parameters, start_date, infection_state));
        }
    }

    // Assign locations round-robin so each location receives an equal share of people.
    int school_ctr = 0, work_ctr = 0, shop_ctr = 0, event_ctr = 0;
    for (auto& person : model.get_persons()) {
        const auto id = person.get_id();
        model.assign_location(id, hospital);
        model.assign_location(id, icu);
        model.assign_location(id, shops[shop_ctr++ % n_shops]);
        model.assign_location(id, events[event_ctr++ % n_events]);
        if (person.get_age() == age_group_5_to_14) {
            model.assign_location(id, schools[school_ctr++ % n_schools]);
        }
        if (person.get_age() == age_group_15_to_34 || person.get_age() == age_group_35_to_59) {
            model.assign_location(id, workplaces[work_ctr++ % n_workplaces]);
        }
    }

    // Set start and end time for the simulation.
    auto t0   = mio::abm::TimePoint(0);
    auto tmax = t0 + mio::abm::days(30);
    auto sim  = mio::abm::Simulation(t0, std::move(model));

    // Both loggers fire at 8AM — one row per school-day morning.
    mio::History<mio::DataWriterToMemory, mio::abm::LogCTCluster, mio::abm::LogSchoolCohort, mio::abm::LogInfectionState>
        historyTimeSeries{
            mio::abm::LogCTCluster{age_group_5_to_14, mio::abm::LocationType::School},
            mio::abm::LogSchoolCohort{},
            mio::abm::LogInfectionState{}};
    // Run the simulation until tmax with the history object.
    sim.advance(tmax, historyTimeSeries);

    // Write results to files.
    auto outpath = mio::create_directories_or_exit(mio::example_results_dir("abm_minimal")) / "ct_cluster.txt";
    std::ofstream outfile(outpath);
    for (const auto& [time, hist] : std::get<0>(historyTimeSeries.get_log())) {
        outfile << time.days();
        for (Eigen::Index i = 0; i < hist.size(); ++i)
            outfile << " " << hist[i];
        outfile << "\n";
    }
    std::cout << "CT cluster log written to " << outpath << std::endl;

    auto cohort_bin_path = mio::create_directories_or_exit(mio::example_results_dir("abm_minimal")) / "cohort.bin";
    auto cohort_txt_path = mio::create_directories_or_exit(mio::example_results_dir("abm_minimal")) / "cohort.txt";
    std::ofstream cohort_bin(cohort_bin_path, std::ios::binary);
    std::ofstream cohort_txt(cohort_txt_path);
    for (const auto& [time, cts] : std::get<1>(historyTimeSeries.get_log())) {
        uint32_t day = static_cast<uint32_t>(std::round(time.days()));
        cohort_bin.write(reinterpret_cast<const char*>(&day), sizeof(day));
        cohort_bin.write(reinterpret_cast<const char*>(cts.data()), cts.size());
        cohort_txt << "day=" << day << ":";
        for (uint8_t ct : cts) {
            cohort_txt << " " << (ct == 255 ? "-" : std::to_string(ct));
        }
        cohort_txt << "\n";
    }
    std::cout << "Cohort log written to " << cohort_bin_path << " and " << cohort_txt_path << std::endl;

    auto inf_path = mio::create_directories_or_exit(mio::example_results_dir("abm_minimal")) / "infection_states.txt";
    std::ofstream inf_file(inf_path);
    inf_file << "t S E I_NS I_Sy I_Sev I_Crit R D\n";
    for (const auto& [time, counts] : std::get<2>(historyTimeSeries.get_log())) {
        inf_file << time.days();
        for (Eigen::Index i = 0; i < counts.size(); ++i) {
            inf_file << " " << counts[i];
        }
        inf_file << "\n";
    }
    std::cout << "Infection state log written to " << inf_path << std::endl;

    return 0;
}
