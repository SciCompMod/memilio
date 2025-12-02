/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding, Sascha Korf
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
#include "abm/result_simulation.h"
#include "abm/household.h"
#include "abm/lockdown_rules.h"
#include "abm/model.h"
#include "abm/time.h"

#include "memilio/compartments/parameter_studies.h"
#include "memilio/data/analyze_result.h"
#include "memilio/io/io.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/base_dir.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/stl_util.h"

#include <string>

constexpr size_t num_age_groups = 4;

/// An ABM setup taken from abm_minimal.cpp.
mio::abm::Model make_model()
{

    const auto age_group_0_to_4   = mio::AgeGroup(0);
    const auto age_group_5_to_14  = mio::AgeGroup(1);
    const auto age_group_15_to_34 = mio::AgeGroup(2);
    const auto age_group_35_to_59 = mio::AgeGroup(3);
    // Create the model with 4 age groups.
    auto model = mio::abm::Model(num_age_groups);

    // Set same infection parameter for all age groups. For example, the incubation period is log normally distributed with parameters 4 and 1.
    model.parameters.get<mio::abm::TimeExposedToNoSymptoms>() = mio::ParameterDistributionLogNormal(4., 1.);

    // Set the age groups that can go to school; here this is AgeGroup(1) (i.e. 5-14)
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age groups that can go to work; here these are AgeGroup(2) and AgeGroup(3) (i.e. 15-34 and 35-59)
    model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({age_group_15_to_34, age_group_35_to_59}, true);

    // Check if the parameters satisfy their constraints.
    model.parameters.check_constraints();

    // There are 10 households for each household group.
    int n_households = 10;

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
    // Maximum contacts limit the number of people that a person can infect while being at this location.
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
    // At every workplace, maximum contacts are 20.
    auto work = model.add_location(mio::abm::LocationType::Work);
    model.get_location(work).get_infection_parameters().set<mio::abm::MaximumContacts>(20);

    // Increase aerosol transmission for all locations
    model.parameters.get<mio::abm::AerosolTransmissionRates>() = 10.0;
    // Increase contact rate for all people between 15 and 34 (i.e. people meet more often in the same location)
    model.get_location(work)
        .get_infection_parameters()
        .get<mio::abm::ContactRates>()[{age_group_15_to_34, age_group_15_to_34}] = 10.0;

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
    // The infection states are chosen randomly with the following distribution
    std::vector<ScalarType> infection_distribution{0.5, 0.3, 0.05, 0.05, 0.05, 0.05, 0.0, 0.0};
    for (auto& person : model.get_persons()) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), infection_distribution));
        auto person_rng = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(person_rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model.parameters, start_date, infection_state));
        }
    }

    // Assign locations to the people
    for (auto& person : model.get_persons()) {
        const auto id = person.get_id();
        //assign shop and event
        model.assign_location(id, event);
        model.assign_location(id, shop);
        //assign hospital and ICU
        model.assign_location(id, hospital);
        model.assign_location(id, icu);
        //assign work/school to people depending on their age
        if (person.get_age() == age_group_5_to_14) {
            model.assign_location(id, school);
        }
        if (person.get_age() == age_group_15_to_34 || person.get_age() == age_group_35_to_59) {
            model.assign_location(id, work);
        }
    }

    // During the lockdown, social events are closed for 90% of people.
    auto t_lockdown = mio::abm::TimePoint(0) + mio::abm::days(10);
    mio::abm::close_social_events(t_lockdown, 0.9, model.parameters);

    return model;
}

int main()
{
    mio::mpi::init();

    mio::set_log_level(mio::LogLevel::warn);

    // Set start and end time for the simulation.
    auto t0   = mio::abm::TimePoint(0);
    auto tmax = t0 + mio::abm::days(5);
    // Set the number of simulations to run in the study
    const size_t num_runs = 3;

    // Create a parameter study.
    // Note that the study for the ABM currently does not make use of the arguments "parameters" or "dt", as we create
    // a new model for each simulation. Hence we set both arguments to 0.
    // This is mostly due to https://github.com/SciCompMod/memilio/issues/1400
    mio::ParameterStudy study(make_model(), t0, tmax, mio::abm::TimeSpan(0), num_runs);

    // Optional: set seeds to get reproducable results
    // study.get_rng().seed({12341234, 53456, 63451, 5232576, 84586, 52345});

    const std::string result_dir = mio::path_join(mio::base_dir(), "example_results");
    if (!mio::create_directory(result_dir)) {
        mio::log_error("Could not create result directory \"{}\".", result_dir);
        return 1;
    }

    auto ensemble_results = study.run(
        [](auto&& model, auto t0_, auto, size_t) {
            auto copy = model;
            copy.reset_rng(mio::thread_local_rng().get_seeds());
            return mio::abm::ResultSimulation(std::move(copy), t0_);
        },
        [result_dir](auto&& sim, auto&& run_idx) {
            auto interpolated_result = mio::interpolate_simulation_result(sim.get_result());
            std::string outpath = mio::path_join(result_dir, "abm_minimal_run_" + std::to_string(run_idx) + ".txt");
            std::ofstream outfile_run(outpath);
            sim.get_result().print_table(outfile_run, {"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4);
            std::cout << "Results written to " << outpath << std::endl;
            auto params = std::vector<mio::abm::Model>{};
            return std::vector{interpolated_result};
        });

    if (ensemble_results.size() > 0) {
        auto ensemble_results_p05 = ensemble_percentile(ensemble_results, 0.05);
        auto ensemble_results_p25 = ensemble_percentile(ensemble_results, 0.25);
        auto ensemble_results_p50 = ensemble_percentile(ensemble_results, 0.50);
        auto ensemble_results_p75 = ensemble_percentile(ensemble_results, 0.75);
        auto ensemble_results_p95 = ensemble_percentile(ensemble_results, 0.95);

        mio::unused(save_result(ensemble_results_p05, {0}, num_age_groups,
                                mio::path_join(result_dir, "Results_" + std::string("p05") + ".h5")));
        mio::unused(save_result(ensemble_results_p25, {0}, num_age_groups,
                                mio::path_join(result_dir, "Results_" + std::string("p25") + ".h5")));
        mio::unused(save_result(ensemble_results_p50, {0}, num_age_groups,
                                mio::path_join(result_dir, "Results_" + std::string("p50") + ".h5")));
        mio::unused(save_result(ensemble_results_p75, {0}, num_age_groups,
                                mio::path_join(result_dir, "Results_" + std::string("p75") + ".h5")));
        mio::unused(save_result(ensemble_results_p95, {0}, num_age_groups,
                                mio::path_join(result_dir, "Results_" + std::string("p95") + ".h5")));
    }

    mio::mpi::finalize();

    return 0;
}
