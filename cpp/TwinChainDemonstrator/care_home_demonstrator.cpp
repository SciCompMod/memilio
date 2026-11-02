/*
* Copyright (C) 2020-2026 MEmilio
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
#include "abm/household.h"
#include "abm/location_type.h"
#include "abm/lockdown_rules.h"
#include "abm/model.h"
#include "abm/common_abm_loggers.h"
#include "abm/parameters.h"
#include "memilio/io/directories.h"

#include <fstream>
#include <vector>
#include <string>

struct CustomLogger : mio::LogAlways {
    using Type = std::vector<
        std::tuple<mio::abm::LocationId, mio::abm::LocationType, mio::abm::PersonId, mio::abm::InfectionState, double>>;
    static Type log(const mio::abm::Simulation<mio::abm::Model>& sim)
    {
        std::vector<std::tuple<mio::abm::LocationId, mio::abm::LocationType, mio::abm::PersonId,
                               mio::abm::InfectionState, double>>
            log_vec{};
        for (auto&& person : sim.get_model().get_persons()) {
            double shed = 0;
            if (person.is_infected(sim.get_time())) {
                shed = person.get_infection().get_viral_shed(sim.get_time());
            }
            log_vec.push_back(std::make_tuple(person.get_location(), person.get_location_type(), person.get_id(),
                                              person.get_infection_state(sim.get_time()), shed));
        }
        return log_vec;
    }
};

void write_custom_log_to_file(const mio::History<mio::DataWriterToMemory, CustomLogger>& history,
                              const std::string& filename, const mio::abm::Simulation<mio::abm::Model>& sim)
{
    auto& log = std::get<0>(history.get_log());
    std::ofstream outfile(filename);
    outfile << "Time,Tau,LocationType,LocationId,PersonId,InfectionState,ViralShed,InfectionHistory" << std::endl;
    for (size_t t = 0; t < log.size(); ++t) {
        for (size_t p = 0; p < log[t].size(); ++p) {
            auto entry = log[t][p];
            if (t > 0 && std::get<1>(entry) == mio::abm::LocationType::Work &&
                std::get<1>(log[t - 1][p]) != mio::abm::LocationType::Work) {
                outfile << t << "," << (17 - (t % 24)) << "," << static_cast<int>(std::get<1>(entry)) << ","
                        << std::get<0>(entry) << "," << std::get<2>(entry) << ","
                        << static_cast<int>(std::get<3>(entry)) << "," << std::get<4>(entry);
                if (std::get<3>(entry) != mio::abm::InfectionState::Susceptible) {
                    const auto& infection_course =
                        sim.get_model().get_person(std::get<2>(entry)).get_infection().get_infection_course();
                    outfile << ",[";
                    for (size_t i = 0; i < infection_course.size(); ++i) {
                        outfile << "(" << infection_course[i].first.hours() << ","
                                << static_cast<int>(infection_course[i].second) << ")";
                        if (i < infection_course.size() - 1) {
                            outfile << ",";
                        }
                    }
                    outfile << "]";
                }
                else {
                    outfile << ",[]";
                }
                outfile << std::endl;
            }
        }
    }
}

int main()
{
    mio::set_log_level(mio::LogLevel::warn);
    // One age group goes to work and the other does not.
    size_t num_age_groups           = 2;
    const auto age_group_go_to_work = mio::AgeGroup(0);
    const auto age_group_other      = mio::AgeGroup(1);

    // Create the model with 2 age groups.
    auto model = mio::abm::Model(num_age_groups);
    // Set same infection parameter for all age groups
    model.parameters.get<mio::abm::TimeExposedToNoSymptoms>()       = mio::ParameterDistributionLogNormal(1.35, 0.25);
    model.parameters.get<mio::abm::SymptomsPerInfectedNoSymptoms>() = 0.6;
    model.parameters.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>() =
        mio::ParameterDistributionLogNormal(0.91, 0.61);
    model.parameters.get<mio::abm::SeverePerInfectedSymptoms>() = 0.0;

    // Set the age group the can go to work
    model.parameters.get<mio::abm::AgeGroupGotoWork>()                       = false;
    model.parameters.get<mio::abm::AgeGroupGotoWork>()[age_group_go_to_work] = true;

    // Check if the parameters satisfy their contraints.
    model.parameters.check_constraints();

    // Only one household with 100 persons, 50 of them go to work and 50 do not.
    size_t n_persons = 100;
    auto home        = model.add_location(mio::abm::LocationType::Home);
    // Add work location
    auto work = model.add_location(mio::abm::LocationType::Work);
    for (size_t i = 0; i < n_persons; ++i) {
        if (i < n_persons / 2) {
            auto person = model.add_person(home, age_group_go_to_work);
            model.get_person(person).set_assigned_location(mio::abm::LocationType::Work, work, model.get_id());
            model.get_person(person).set_assigned_location(mio::abm::LocationType::Home, home, model.get_id());
        }
        else {
            auto person = model.add_person(home, age_group_other);
            model.get_person(person).set_assigned_location(mio::abm::LocationType::Home, home, model.get_id());
        }
    }

    // Consider only contact transmission
    model.parameters.get<mio::abm::AerosolTransmissionRates>() = 0;
    // Set contact rate at home
    model.get_location(home).get_infection_parameters().get<mio::abm::ContactRates>().get_baseline().setConstant(3.0);
    // Set infection rate from viral shed
    model.parameters.get<mio::abm::InfectionRateFromViralShed>() = 0.5;

    // Assign infection state to each person.
    // The infection states are chosen randomly with the following distribution
    auto start_date = mio::abm::TimePoint(0);
    std::vector<ScalarType> infection_distribution{0.9, 0.05, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (auto& person : model.get_persons()) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), infection_distribution));
        auto rng = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model.parameters, start_date, infection_state));
        }
    }

    // Run the simulation for 10 days.
    auto tmax = start_date + mio::abm::days(10);
    auto sim  = mio::abm::Simulation(start_date, std::move(model));

    // Create a history object to store the time series of the infection states.
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)};
    mio::History<mio::DataWriterToMemory, CustomLogger> historyCustom;

    // Run the simulation until tmax with the history object.
    sim.advance(tmax, historyTimeSeries, historyCustom);

    write_custom_log_to_file(
        historyCustom,
        mio::create_directories_or_exit(mio::example_results_dir("care_home_demonstrator")) / "output.csv", sim);

    // Write results to a file. Also print the filepath to make it easier to find
    auto outpath = mio::create_directories_or_exit(mio::example_results_dir("care_home_demonstrator")) / "history.txt";
    std::ofstream outfile(outpath);
    std::get<0>(historyTimeSeries.get_log())
        .print_table(outfile, {"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4);
    std::cout << "Results written to " << outpath << std::endl;

    return 0;
}
