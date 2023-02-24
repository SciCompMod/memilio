/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include "ode_seair/model.h"
#include "ode_seair/infection_state.h"
#include "ode_seair/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/math/adapt_rk.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/time_series_to_file.h"
#include <fstream>



int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 100;
    double dt   = 0.2;

    mio::log_info("Simulating SEAIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::oseair::Model model;
    const double N = 327167434;// total population of the US

    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}] = 0.9977558755803503;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}]   = 0.0003451395725394549;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}]   = 0.00037846880968213874;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}]  = (337072.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}] = (17448.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Perished)}]   = (9619.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::ObjectiveFunction)}]   = 0.0;




    model.check_constraints();
    // print_seir_params(model);
    auto integrator = std::make_shared<mio::RKIntegratorCore>();
    integrator->set_dt_max(dt);
    integrator->set_abs_tolerance(1e-6);
    integrator->set_rel_tolerance(1e-6);

    auto seair = simulate(t0, tmax, dt, model, integrator);
    const std::string file_name = "seair.txt";
    std::cout << "Writing output to " << file_name << std::endl;
    mio::time_series_to_file(seair, file_name);

}
