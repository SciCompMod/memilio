/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn
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

#include "ode_sir/model.h"
#include "ode_sir/infection_state.h"
#include "ode_sir/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0.;
    double tmax = 50.;
    double dt   = 0.1002004008016032;

    double total_population = 1061000;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::osir::Model model(1);

    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Infected}]  = 1000;
    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Recovered}] = 1000;
    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Susceptible}] =
        total_population -
        model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Recovered}];
    model.parameters.set<mio::osir::TimeInfected>(2);
    model.parameters.set<mio::osir::TransmissionProbabilityOnContact>(1);
    model.parameters.get<mio::osir::ContactPatterns>().get_cont_freq_mat()[0].get_baseline().setConstant(2.7);
    model.parameters.get<mio::osir::ContactPatterns>().get_cont_freq_mat()[0].add_damping(0.6, mio::SimulationTime(12.5));

    auto integrator = std::make_shared<mio::EulerIntegratorCore>();

    model.check_constraints();

    auto sir = simulate(t0, tmax, dt, model, integrator);

    bool print_to_terminal = true;

    if (print_to_terminal) {
        std::vector<std::string> vars = {"S", "I", "R"};
        printf("\n # t");
        for (size_t k = 0; k < (size_t)mio::osir::InfectionState::Count; k++) {
            printf(" %s", vars[k].c_str());
        }

        auto num_points = static_cast<size_t>(sir.get_num_time_points());
        for (size_t i = 0; i < num_points; i++) {
            printf("\n%.14f ", sir.get_time(i));
            Eigen::VectorXd res_j = sir.get_value(i);
            for (size_t j = 0; j < (size_t)mio::osir::InfectionState::Count; j++) {
                printf(" %.14f", res_j[j]);
            }
        }

        Eigen::VectorXd res_j = sir.get_last_value();
        printf("number total: %f", res_j[0] + res_j[1] + res_j[2]);
    }
}
