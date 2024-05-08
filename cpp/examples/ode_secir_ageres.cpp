/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#include "ode_secir/model.h"
#include "ode_secir/parameters_io.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/compartments/simulation.h"

int main()
{

    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    mio::log_info("Simulating SECIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    mio::osecir::Model model(6);
    auto nb_groups = model.parameters.get_num_groups();
    double fact    = 1.0 / (double)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<mio::osecir::StartDay>(60);
    params.set<mio::osecir::Seasonality>(0.2);
    params.get<mio::osecir::TestAndTraceCapacity>() = 35;

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::TimeExposed>()[i]            = 3.2;
        params.get<mio::osecir::TimeInfectedNoSymptoms>()[i] = 2.;
        params.get<mio::osecir::TimeInfectedSymptoms>()[i]   = 5.8;
        params.get<mio::osecir::TimeInfectedSevere>()[i]     = 9.5;
        params.get<mio::osecir::TimeInfectedCritical>()[i]   = 7.1;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]                     = fact * nb_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}]          = fact * nb_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]            = fact * nb_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]              = fact * nb_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]            = fact * nb_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]                   = fact * nb_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]                        = fact * nb_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * nb_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact>()[i]  = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]    = 0.7;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]    = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]    = 0.25;
        params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>()[i] = 0.45;
        params.get<mio::osecir::SeverePerInfectedSymptoms>()[i]         = 0.2;
        params.get<mio::osecir::CriticalPerSevere>()[i]                 = 0.45;
        params.get<mio::osecir::DeathsPerCritical>()[i]                 = 0.3;
    }

    model.apply_constraints();

    // tests
    const auto num_age_groups               = 6;
    const auto TEST_DATA_DIR                = "/localdata1/code_2024/memilio/data/";
    const auto results_dir                  = "/localdata1/code_2024/memilio/test";
    std::vector<mio::osecir::Model> models1 = {model};
    auto status_export                      = mio::osecir::export_input_data_county_timeseries(
        models1, results_dir, {1001}, {2020, 12, 01}, std::vector<double>(size_t(num_age_groups), 1.0), 1.0, 1,
        mio::path_join(TEST_DATA_DIR, "pydata/Germany", "county_divi_ma7.json"),
        mio::path_join(TEST_DATA_DIR, "pydata/Germany", "cases_all_county_age_ma7.json"),
        mio::path_join(TEST_DATA_DIR, "pydata/Germany", "county_current_population.json"));

    std::vector<mio::osecir::Model> models2 = {model};
    auto read_result1                       = mio::osecir::read_input_data_county(models2, {2020, 12, 01}, {1001},
                                                            std::vector<double>(size_t(num_age_groups), 1.0), 1.0,
                                                            TEST_DATA_DIR, 1, false);

    mio::TimeSeries<double> secir = simulate(0.0, 1.0, 1.0, models2[0]);

    std::cout << secir.get_value(0) << std::endl;

    mio::unused(read_result1);
}
