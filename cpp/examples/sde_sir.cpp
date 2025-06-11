/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Nils Wassmuth, Rene Schmieding, Martin J. Kuehn
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

#include "memilio/compartments/stochastic_simulation.h"
#include "sde_sir/model.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0.;
    double tmax = 5.;
    double dt   = 0.1;

    double total_population = 10000;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    // using Model = ;
    mio::ssir::Model model;

    model.populations[{mio::Index<mio::ssir::InfectionState>(mio::ssir::InfectionState::Infected)}]  = 100;
    model.populations[{mio::Index<mio::ssir::InfectionState>(mio::ssir::InfectionState::Recovered)}] = 1000;
    model.populations[{mio::Index<mio::ssir::InfectionState>(mio::ssir::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::ssir::InfectionState>(mio::ssir::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::ssir::InfectionState>(mio::ssir::InfectionState::Recovered)}];
    model.parameters.set<mio::ssir::TimeInfected>(10);
    model.parameters.set<mio::ssir::TransmissionProbabilityOnContact>(1);
    model.parameters.get<mio::ssir::ContactPatterns>().get_baseline()(0, 0) = 2.7;
    model.parameters.get<mio::ssir::ContactPatterns>().add_damping(0.6, mio::SimulationTime(12.5));

    model.check_constraints();

    auto sir = mio::simulate_stochastic(t0, tmax, dt, model);

    sir.print_table();
}
