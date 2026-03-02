/*
* Copyright (C) 2020-2026 MEmilio
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

    ScalarType t0   = 0;
    ScalarType tmax = 50;
    ScalarType dt   = 0.1;

    mio::log_info("Simulating SECIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    ScalarType cont_freq = 10; // see Polymod study

    ScalarType nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
               nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model<ScalarType> model(3);
    auto nb_groups  = model.parameters.get_num_groups();
    ScalarType fact = 1.0 / (ScalarType)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<mio::osecir::StartDay<ScalarType>>(60);
    params.set<mio::osecir::Seasonality<ScalarType>>(0.2);
    params.get<mio::osecir::TestAndTraceCapacity<ScalarType>>() = 35;

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::TimeExposed<ScalarType>>()[i]            = 3.2;
        params.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[i] = 2.;
        params.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[i]   = 5.8;
        params.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[i]     = 9.5;
        params.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[i]   = 7.1;

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

        params.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[i]  = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[i]    = 0.7;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i]    = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[i]    = 0.25;
        params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<ScalarType>>()[i] = 0.45;
        params.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[i]         = 0.2;
        params.get<mio::osecir::CriticalPerSevere<ScalarType>>()[i]                 = 0.25;
        params.get<mio::osecir::DeathsPerCritical<ScalarType>>()[i]                 = 0.3;
    }
    // The function apply_constraints() ensures that all parameters are within their defined bounds.
    // Note that negative values are set to zero instead of stopping the simulation.
    model.apply_constraints();

    mio::ContactMatrixGroup<ScalarType>& contact_matrix = params.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0]                                   = mio::ContactMatrix<ScalarType>(
        Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                               mio::SimulationTime<ScalarType>(30.));

    mio::TimeSeries<ScalarType> secir = mio::simulate<ScalarType, mio::osecir::Model<ScalarType>>(t0, tmax, dt, model);
    bool print_to_terminal            = true;

    if (print_to_terminal) {

        std::vector<std::string> vars = {"S", "E", "C", "C_confirmed", "I", "I_confirmed", "H", "U", "R", "D"};
        printf("Number of time points :%d\n", static_cast<int>(secir.get_num_time_points()));
        printf("People in\n");

        for (size_t k = 0; k < (size_t)mio::osecir::InfectionState::Count; k++) {
            ScalarType dummy = 0;

            for (size_t i = 0; i < (size_t)params.get_num_groups(); i++) {
                printf("\t %s[%d]: %.0f", vars[k].c_str(), (int)i,
                       secir.get_last_value()[k + (size_t)mio::osecir::InfectionState::Count * (int)i]);
                dummy += secir.get_last_value()[k + (size_t)mio::osecir::InfectionState::Count * (int)i];
            }

            printf("\t %s_otal: %.0f\n", vars[k].c_str(), dummy);
        }
    }
}
