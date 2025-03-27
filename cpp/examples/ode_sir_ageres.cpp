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
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/compartments/simulation.h"
#include "ode_sir/model.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    double cont_freq = 10; // see Polymod study

    double nb_total_t0 = 10000, nb_inf_t0 = 50, nb_rec_t0 = 10;

    const size_t num_groups = 3;
    mio::osir::Model model(num_groups);

    double fact = 1.0 / num_groups;

    auto& params = model.parameters;

    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(num_groups); i++) {
        model.populations[{i, mio::osir::InfectionState::Infected}]  = fact * nb_inf_t0;
        model.populations[{i, mio::osir::InfectionState::Recovered}] = fact * nb_rec_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osir::InfectionState::Susceptible},
                                                                         fact * nb_total_t0);

        model.parameters.get<mio::osir::TimeInfected<double>>()[i]                     = 2.0;
        model.parameters.get<mio::osir::TransmissionProbabilityOnContact<double>>()[i] = 0.3;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osir::ContactPatterns<double>>();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, fact * cont_freq));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant(num_groups, num_groups, 0.7), mio::SimulationTime(30.));

    model.apply_constraints();

    auto sir = simulate(t0, tmax, dt, model);

    std::vector<std::string> vars = {"S", "I", "R"};
    printf("Number of time points :%d\n", static_cast<int>(sir.get_num_time_points()));
    printf("People in\n");

    for (size_t k = 0; k < (size_t)mio::osir::InfectionState::Count; k++) {
        double dummy = 0;

        for (size_t i = 0; i < (size_t)params.get_num_groups(); i++) {
            printf("\t %s[%d]: %.0f", vars[k].c_str(), (int)i,
                   sir.get_last_value()[k + (size_t)mio::osir::InfectionState::Count * (int)i]);
            dummy += sir.get_last_value()[k + (size_t)mio::osir::InfectionState::Count * (int)i];
        }

        printf("\t %s_Total: %.0f\n", vars[k].c_str(), dummy);
    }
}