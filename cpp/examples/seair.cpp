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
#include "ad/ad.hpp"



int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    using FP=double;



    FP  t0   = 0;
    FP  tmax = 100;
    FP  dt   = 0.2;

    mio::log_info("Simulating SEAIR; t={} ... {} with dt = {}.", ad::value(t0), ad::value(tmax), ad::value(dt));

    mio::oseair::Model<FP> model;
    const double N = 327167434;// total population of the United States

    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}] = 0.9977558755803503;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}]   = 0.0003451395725394549;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}]   = 0.00037846880968213874;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}]  = (337072.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}] = (17448.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Perished)}]   = (9619.0 / N);
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::ObjectiveFunction)}]   = 0.0;




    model.check_constraints();




    auto seair = mio::simulate<mio::oseair::Model<FP>,FP>(t0, tmax, dt, model);
//    const std::string file_name = "seair.txt";
//    std::cout << "Writing output to " << file_name << std::endl;
//    mio::time_series_to_file(seair, file_name);

}
