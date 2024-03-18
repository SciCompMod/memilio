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

#include "sde_sir/model.h"
#include "sde_sir/infection_state.h"
#include "sde_sir/parameters.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/random_number_generator.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0.;
    double tmax = 5.;
    double dt   = 0.1;

    double total_population = 10000;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

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

    auto integrator = std::make_shared<mio::EulerIntegratorCore>();

    model.check_constraints();


    //auto sir = simulate(t0, tmax, dt, model, integrator);
    /*auto integrator = std::make_shared<mio::EulerIntegratorCore>();
    
    sim.set_integrator(integrator);
      */ 
    //mio::Simulation<model> sim(model, t0, dt, integrator);
    //sim.advance(tmax); 
    mio::FlowSimulation<mio::ssir::Model> sim(model, t0, dt);
    sim.set_integrator(integrator);
    //auto sir = simulate_stoch(t0, tmax, dt, model, integrator);
    sim.advance_stochastic(tmax);
    auto sir = sim.get_result();
    bool print_to_terminal = true;

    sir.print_table();
}
