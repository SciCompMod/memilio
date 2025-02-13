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

#include "memilio/utils/logging.h"
#include "sde_sirs/model.h"
#include "sde_sirs/simulation.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0.;
    double tmax = 5.;
    double dt   = 0.001;

    double total_population = 10000;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::ssirs::Model model;

    model.populations[{mio::Index<mio::ssirs::InfectionState>(mio::ssirs::InfectionState::Infected)}]  = 100;
    model.populations[{mio::Index<mio::ssirs::InfectionState>(mio::ssirs::InfectionState::Recovered)}] = 1000;
    model.populations[{mio::Index<mio::ssirs::InfectionState>(mio::ssirs::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::ssirs::InfectionState>(mio::ssirs::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::ssirs::InfectionState>(mio::ssirs::InfectionState::Recovered)}];
    model.parameters.set<mio::ssirs::TimeInfected>(10);
    model.parameters.set<mio::ssirs::TimeImmune>(100);
    model.parameters.set<mio::ssirs::TransmissionProbabilityOnContact>(1);
    model.parameters.get<mio::ssirs::ContactPatterns>().get_baseline()(0, 0) = 20.7;
    model.parameters.get<mio::ssirs::ContactPatterns>().add_damping(0.6, mio::SimulationTime(12.5));

    model.check_constraints();

    auto ssirs = mio::ssirs::simulate(t0, tmax, dt, model);

    ssirs.print_table({"Susceptible", "Infected", "Recovered"});
}
