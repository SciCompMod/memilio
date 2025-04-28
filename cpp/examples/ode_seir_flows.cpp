/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Martin J. Kuehn, Rene Schmieding, Henrik Zunker
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
#include "memilio/compartments/flow_simulation.h"
#include "memilio/utils/logging.h"

#include <vector>

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 1;
    double dt   = 0.001;

    mio::log_info("Simulating SEIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::oseir::Model<double> model(1);

    constexpr double total_population                                            = 10000;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]   = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]  = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}] = 100;
    model.populations.set_difference_from_group_total<mio::AgeGroup>(
        {mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}, total_population);

    model.parameters.set<mio::oseir::TimeExposed<double>>(5.2);
    model.parameters.set<mio::oseir::TimeInfected<double>>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<double>>(0.04);

    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::oseir::ContactPatterns<double>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(10);

    model.check_constraints();

    auto seir = simulate_flows(t0, tmax, dt, model);

    printf("Compartments: \n");
    seir[0].print_table({"S", "E", "I", "R"});

    printf("Flows: \n");
    seir[1].print_table({"S->E", "E->I", "I->R"});

    printf("\n number total: %f\n", seir[0].get_last_value()[0] + seir[0].get_last_value()[1] +
                                        seir[0].get_last_value()[2] + seir[0].get_last_value()[3]);
}
