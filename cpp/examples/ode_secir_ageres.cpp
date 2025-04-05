/* 
* Copyright (C) 2020-2025 MEmilio
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

    double cont_freq = 10; // see Polymod study

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    mio::osecir::Model<double> model(3);
    auto nb_groups = model.parameters.get_num_groups();
    double fact    = 1.0 / (double)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<mio::osecir::StartDay>(60);
    params.set<mio::osecir::Seasonality<double>>(0.2);
    params.get<mio::osecir::TestAndTraceCapacity<double>>() = 35;

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::TimeExposed<double>>()[i]            = 3.2;
        params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i] = 2.;
        params.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]   = 5.8;
        params.get<mio::osecir::TimeInfectedSevere<double>>()[i]     = 9.5;
        params.get<mio::osecir::TimeInfectedCritical<double>>()[i]   = 7.1;

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

        params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i]  = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[i]    = 0.7;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i]    = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i]    = 0.25;
        params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>()[i] = 0.45;
        params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i]         = 0.2;
        params.get<mio::osecir::CriticalPerSevere<double>>()[i]                 = 0.25;
        params.get<mio::osecir::DeathsPerCritical<double>>()[i]                 = 0.3;
    }

    model.apply_constraints();

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                               mio::SimulationTime(30.));

    mio::TimeSeries<double> secir = mio::simulate<double, mio::osecir::Model<double>>(t0, tmax, dt, model);
    bool print_to_terminal        = true;

    if (print_to_terminal) {

        std::vector<std::string> vars = {"S", "E", "C", "C_confirmed", "I", "I_confirmed", "H", "U", "R", "D"};
        printf("Number of time points :%d\n", static_cast<int>(secir.get_num_time_points()));
        printf("People in\n");

        for (size_t k = 0; k < (size_t)mio::osecir::InfectionState::Count; k++) {
            double dummy = 0;

            for (size_t i = 0; i < (size_t)params.get_num_groups(); i++) {
                printf("\t %s[%d]: %.0f", vars[k].c_str(), (int)i,
                       secir.get_last_value()[k + (size_t)mio::osecir::InfectionState::Count * (int)i]);
                dummy += secir.get_last_value()[k + (size_t)mio::osecir::InfectionState::Count * (int)i];
            }

            printf("\t %s_otal: %.0f\n", vars[k].c_str(), dummy);
        }
    }
}
