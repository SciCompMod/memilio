/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Wadim Koslow
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
#include "memilio/io/result_io.h"

#include <iostream>

int main()
{

    const auto t0   = 0.;
    const auto tmax = 10.;
    const auto dt   = 1.; //time step of mobility, not integration

    double cont_freq = 10; // see Polymod study

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model<double> model(1);
    mio::AgeGroup nb_groups = model.parameters.get_num_groups();
    double fact             = 1.0 / (double)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<mio::osecir::ICUCapacity<double>>(std::numeric_limits<double>::max());
    params.set<mio::osecir::StartDay<double>>(0);
    params.set<mio::osecir::Seasonality<double>>(0);

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::TimeExposed<double>>()[i]            = 3.2;
        params.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[i] = 2.0;
        params.get<mio::osecir::TimeInfectedSymptoms<double>>()[i]   = 6.;
        params.get<mio::osecir::TimeInfectedSevere<double>>()[i]     = 12;
        params.get<mio::osecir::TimeInfectedCritical<double>>()[i]   = 8;

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

        params.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[i] = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[i]   = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()[i]   = 0.25;
        params.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[i]        = 0.2;
        params.get<mio::osecir::CriticalPerSevere<double>>()[i]                = 0.25;
        params.get<mio::osecir::DeathsPerCritical<double>>()[i]                = 0.3;
    }

    params.apply_constraints();

    mio::ContactMatrixGroup<double>& contact_matrix = params.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] =
        mio::ContactMatrix<double>(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    contact_matrix.add_damping(0.3, mio::SimulationTime<double>(30.));

    auto result_from_sim                                  = mio::osecir::simulate<double>(t0, tmax, dt, model);
    std::vector<mio::TimeSeries<double>> results_from_sim = {result_from_sim, result_from_sim};
    std::vector<int> ids                                  = {1, 2};

    auto save_result_status = mio::save_result(results_from_sim, ids, (int)(size_t)nb_groups, "test_result.h5");
}
