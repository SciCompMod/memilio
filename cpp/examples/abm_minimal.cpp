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
#include "abm/model.h"
#include "abm/common_abm_loggers.h"

#include <fstream>
#include <iostream>
#include <cstdlib>
std::pair<double, double> get_my_and_sigma(std::pair<double, double> mean_and_std)
{
    auto mean    = mean_and_std.first;
    auto stddev  = mean_and_std.second;
    double my    = log(mean * mean / sqrt(mean * mean + stddev * stddev));
    double sigma = sqrt(log(1 + stddev * stddev / (mean * mean)));
    return {my, sigma};
}
void set_parameters(mio::abm::Parameters& params, mio::AgeGroup age_group_0_to_4)
{
    auto incubation_period_my_sigma          = get_my_and_sigma({4.5, 1.5});
    params.get<mio::abm::TimeExposedToNoSymptoms>() = mio::ParameterDistributionLogNormal(incubation_period_my_sigma.first, incubation_period_my_sigma.second);

    auto InfectedNoSymptoms_to_symptoms_my_sigma             = get_my_and_sigma({1.1, 0.9});
    params.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>() = mio::ParameterDistributionLogNormal(InfectedNoSymptoms_to_symptoms_my_sigma.first,
                                                                InfectedNoSymptoms_to_symptoms_my_sigma.second);

    auto TimeInfectedNoSymptomsToRecovered_my_sigma           = get_my_and_sigma({8.0, 2.0});
    params.get<mio::abm::TimeInfectedNoSymptomsToRecovered>() = mio::ParameterDistributionLogNormal(TimeInfectedNoSymptomsToRecovered_my_sigma.first,
                                                                 TimeInfectedNoSymptomsToRecovered_my_sigma.second);

    auto TimeInfectedSymptomsToSevere_my_sigma           = get_my_and_sigma({6.6, 4.9});
    params.get<mio::abm::TimeInfectedSymptomsToSevere>() = mio::ParameterDistributionLogNormal(TimeInfectedSymptomsToSevere_my_sigma.first,
                                                            TimeInfectedSymptomsToSevere_my_sigma.second);

    auto TimeInfectedSymptomsToRecovered_my_sigma           = get_my_and_sigma({8.0, 2.0});
    params.get<mio::abm::TimeInfectedSymptomsToRecovered>() = mio::ParameterDistributionLogNormal(TimeInfectedSymptomsToRecovered_my_sigma.first,
                                                               TimeInfectedSymptomsToRecovered_my_sigma.second);

    auto TimeInfectedSevereToCritical_my_sigma           = get_my_and_sigma({1.5, 2.0});
    params.get<mio::abm::TimeInfectedSevereToCritical>() = mio::ParameterDistributionLogNormal(TimeInfectedSevereToCritical_my_sigma.first,
                                                            TimeInfectedSevereToCritical_my_sigma.second);

    auto TimeInfectedSevereToRecovered_my_sigma           = get_my_and_sigma({18.1, 6.3});
    params.get<mio::abm::TimeInfectedSevereToRecovered>() = mio::ParameterDistributionLogNormal(TimeInfectedSevereToRecovered_my_sigma.first,
                                                             TimeInfectedSevereToRecovered_my_sigma.second);

    auto TimeInfectedCriticalToDead_my_sigma           = get_my_and_sigma({10.7, 4.8});
    params.get<mio::abm::TimeInfectedCriticalToDead>() = mio::ParameterDistributionLogNormal(TimeInfectedCriticalToDead_my_sigma.first,
                                                          TimeInfectedCriticalToDead_my_sigma.second);

    auto TimeInfectedCriticalToRecovered_my_sigma           = get_my_and_sigma({18.1, 6.3});
    params.get<mio::abm::TimeInfectedCriticalToRecovered>() = mio::ParameterDistributionLogNormal(TimeInfectedCriticalToRecovered_my_sigma.first,
                                                               TimeInfectedCriticalToRecovered_my_sigma.second);

    // Set percentage parameters
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.50;

    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.02;

    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.1;

    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.12;
}
int main(int argc, char* argv[])
{
    // Parse command-line arguments
    // Usage: ./abm_minimal_example [infection_rate] [num_household_members]
    double infection_rate = 10.0; // default value
    int num_household_members = 2; // default value
    
    if (argc > 1) {
        infection_rate = std::atof(argv[1]);
    }
    if (argc > 2) {
        num_household_members = std::atoi(argv[2]);
    }
    
    std::cout << "Running simulation with:" << std::endl;
    std::cout << "  InfectionRateFromViralShed: " << infection_rate << std::endl;
    std::cout << "  Number of household members: " << num_household_members << std::endl;
    
    // This is a minimal example with children and adults < 60 year old.
    // We divided them into 4 different age groups, which are defined as follows:
    mio::set_log_level(mio::LogLevel::warn);
    size_t num_age_groups         = 1;
    const auto age_group_0_to_4   = mio::AgeGroup(0);

    // Create the model with 4 age groups.
    auto model = mio::abm::Model(num_age_groups);
    auto start_date = mio::abm::TimePoint(0);

    for (auto& loc : model.get_locations()) {
        if (loc.get_type()== mio::abm::LocationType::Home) {
             loc.get_infection_parameters().get<mio::abm::ContactRates>()[{age_group_0_to_4, age_group_0_to_4}] = 0.8806*1.6;
        }
    }

    // Set the age group the can go to school is AgeGroup(1) (i.e. 5-14)
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    // Set the age group the can go to work is AgeGroup(2) and AgeGroup(3) (i.e. 15-34 and 35-59)
    model.parameters.get<mio::abm::AgeGroupGotoWork>() = false;

    // Check if the parameters satisfy their contraints.
    model.parameters.check_constraints();

    // There are 10 households for each household group.
    int n_households = 1;
    auto person = mio::abm::HouseholdMember(num_age_groups); // A parent is 50/50% 15-34 or 35-59.
    person.set_age_weight(age_group_0_to_4, 1);

    // Two-person household with one parent and one child.
    auto household  = mio::abm::Household();
    auto hh_group   = mio::abm::HouseholdGroup();
    household.add_members(person, num_household_members);
    hh_group.add_households(household, n_households);
    add_household_group_to_model(model,hh_group);


    // Increase aerosol transmission for all locations
    model.parameters.get<mio::abm::AerosolTransmissionRates>() = 00.0;
    set_parameters(model.parameters, age_group_0_to_4);

    auto first_person = model.get_persons().begin();
    auto rng = mio::abm::PersonalRandomNumberGenerator(*first_person);
    mio::abm::InfectionState infection_state = mio::abm::InfectionState::Exposed;
    auto infection = mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, first_person->get_age(),
                                                         model.parameters, start_date- mio::abm::days(0), infection_state);
    first_person->add_new_infection(std::move(infection));
   
    model.parameters.get<mio::abm::InfectionRateFromViralShed>()[{mio::abm::VirusVariant::Wildtype}] = infection_rate;


    // During the lockdown, social events are closed for 90% of people.
    // auto t_lockdown = mio::abm::TimePoint(0) + mio::abm::days(10);
    // mio::abm::close_social_events(t_lockdown, 0.9, model.parameters);

    // Set start and end time for the simulation.
    auto t0   = mio::abm::TimePoint(0);
    auto tmax = t0 + mio::abm::days(10);
    auto sim  = mio::abm::Simulation(t0, std::move(model));

    // Create a history object to store the time series of the infection states.
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)};

    // Run the simulation until tmax with the history object.
    sim.advance(tmax, historyTimeSeries);

    // The results are written into the file "abm_minimal.txt" as a table with 9 columns.
    // The first column is Time. The other columns correspond to the number of people with a certain infection state at this Time:
    // Time = Time in days, S = Susceptible, E = Exposed, I_NS = InfectedNoSymptoms, I_Sy = InfectedSymptoms, I_Sev = InfectedSevere,
    // I_Crit = InfectedCritical, R = Recovered, D = Dead
    std::ofstream outfile("abm_minimal.txt");
    std::get<0>(historyTimeSeries.get_log())
        .print_table(outfile, {"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4);
    std::cout << "Results written to abm_minimal.txt" << std::endl;

    return 0;
}
