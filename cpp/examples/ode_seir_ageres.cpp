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
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 50;
    double dt   = 0.001;

    mio::log_info("Simulating SEIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    double cont_freq = 10;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_rec_t0 = 10;
    const size_t num_groups = 3;

    mio::oseir::Model<double> model(num_groups);
    double fact = 1.0 / num_groups;

    auto& params = model.parameters;

    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(num_groups); i++) {
        model.populations[{i, mio::oseir::InfectionState::Exposed}]   = fact * nb_exp_t0;
        model.populations[{i, mio::oseir::InfectionState::Infected}]  = fact * nb_inf_t0;
        model.populations[{i, mio::oseir::InfectionState::Recovered}] = fact * nb_rec_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::oseir::InfectionState::Susceptible},
                                                                         fact * nb_total_t0);

        model.parameters.get<mio::oseir::TimeExposed<double>>()[i]                      = 5.2;
        model.parameters.get<mio::oseir::TimeInfected<double>>()[i]                     = 6;
        model.parameters.get<mio::oseir::TransmissionProbabilityOnContact<double>>()[i] = 0.04;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::oseir::ContactPatterns<double>>();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, fact * cont_freq));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant(num_groups, num_groups, 0.7), mio::SimulationTime(30.));

    model.apply_constraints();

    auto seir = simulate(t0, tmax, dt, model);

    std::vector<std::string> vars = {"S", "E", "I", "R"};
    printf("Number of time points :%d\n", static_cast<int>(seir.get_num_time_points()));
    printf("People in\n");

    for (size_t k = 0; k < (size_t)mio::oseir::InfectionState::Count; k++) {
        double dummy = 0;

        for (size_t i = 0; i < (size_t)params.get_num_groups(); i++) {
            printf("\t %s[%d]: %.0f", vars[k].c_str(), (int)i,
                   seir.get_last_value()[k + (size_t)mio::oseir::InfectionState::Count * (int)i]);
            dummy += seir.get_last_value()[k + (size_t)mio::oseir::InfectionState::Count * (int)i];
        }

        printf("\t %s_Total: %.0f\n", vars[k].c_str(), dummy);
    }
}
