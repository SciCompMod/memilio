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
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/math/euler.h"
#include "memilio/utils/time_series.h"

#include "memilio/utils/time_series.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0;
    ScalarType tmax = 0.2;
    ScalarType dt   = 0.1;

    mio::log_info("Simulating ODE SEIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    int number_agegroups = 6;
    mio::oseir::Model<ScalarType> model(number_agegroups);

    ScalarType total_population = 286922;
    for (int i = 0; i < number_agegroups; i++) {
        model.populations[{mio::AgeGroup(i), mio::oseir::InfectionState::Exposed}]   = 10;
        model.populations[{mio::AgeGroup(i), mio::oseir::InfectionState::Infected}]  = 0;
        model.populations[{mio::AgeGroup(i), mio::oseir::InfectionState::Recovered}] = 0;
        model.populations[{mio::AgeGroup(i), mio::oseir::InfectionState::Susceptible}] =
            total_population - model.populations[{mio::AgeGroup(i), mio::oseir::InfectionState::Exposed}] -
            model.populations[{mio::AgeGroup(i), mio::oseir::InfectionState::Infected}] -
            model.populations[{mio::AgeGroup(i), mio::oseir::InfectionState::Recovered}];
    }

    model.parameters.set<mio::oseir::TimeExposed<ScalarType>>(3.);
    model.parameters.set<mio::oseir::TimeInfected<ScalarType>>(5.);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(1.);

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(7.95);
    // contact_matrix[0].add_damping(0.7, mio::SimulationTime(30.));

    model.check_constraints();
    std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator = std::make_shared<mio::EulerIntegratorCore<>>();

    auto seir = simulate(t0, tmax, dt, model, integrator);

    auto reproduction_numbers = model.get_reproduction_numbers(seir);
    std::cout << "\nbasis reproduction number: " << reproduction_numbers[0] << "\n";

    seir.print_table({"S", "E", "I", "R"});
    // std::cout << "\nnumber total: " << seir.get_last_value().sum() << "\n";
}
