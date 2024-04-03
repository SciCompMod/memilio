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

#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "ode_sir/infection_state.h"
#include "ode_sir/model.h"
#include "ode_sir/parameters.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0.;
    double tmax = 50.;
    double dt   = 0.1;

    double total_population = 1061000;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::osir::Model<double> model;

    model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Infected)}]  = 1000;
    model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Recovered)}] = 1000;
    model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Recovered)}];
    model.parameters.set<mio::osir::TimeInfected<double>>(2);
    model.parameters.set<mio::osir::TransmissionProbabilityOnContact<double>>(0.04);
    model.parameters.get<mio::osir::ContactPatterns>().get_baseline()(0, 0) = 2.7;
    model.parameters.get<mio::osir::ContactPatterns>().add_damping(0.6, mio::SimulationTime(12.5));

    auto integrator = std::make_shared<mio::EulerIntegratorCore<double>>();

    model.check_constraints();

    auto sir = mio::simulate<double, mio::osir::Model<double>>(t0, tmax, dt, model, integrator);

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
