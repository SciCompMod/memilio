/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/data/analyze_result.h"
#include "memilio/compartments/simulation.h"
#include "memilio/config.h"
#include "memilio/math/euler.h"
#include "memilio/math/integrator.h"
#include "memilio/utils/logging.h"
#include "ode_sir/infection_state.h"
#include "ode_sir/model.h"
#include "ode_sir/parameters.h"
#include <memory>

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0.;
    ScalarType tmax = 50.;
    ScalarType dt   = 0.1;

    ScalarType total_population = 1061000;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::osir::Model<ScalarType> model(1);

    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Infected}]  = 1000;
    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Recovered}] = 1000;
    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Susceptible}] =
        total_population - model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Infected}] -
        model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Recovered}];
    model.parameters.set<mio::osir::TimeInfected<ScalarType>>(2);
    model.parameters.set<mio::osir::TransmissionProbabilityOnContact<ScalarType>>(0.5);

    mio::ContactMatrixGroup& contact_matrix =
        model.parameters.get<mio::osir::ContactPatterns<ScalarType>>().get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(2.7);
    contact_matrix[0].add_damping(0.6, mio::SimulationTime(12.5));
    model.check_constraints();

    std::shared_ptr<mio::OdeIntegratorCore<ScalarType>> integrator =
        std::make_shared<mio::EulerIntegratorCore<ScalarType>>();
    auto sir = simulate(t0, tmax, dt, model, integrator);

    // interpolate results
    auto interpolated_results = mio::interpolate_simulation_result(sir);

    interpolated_results.print_table({"S", "I", "R"});
    std::cout << "\nPopulation total: " << sir.get_last_value().sum() << "\n";
}
