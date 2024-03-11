/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Anna Wendler, Lena Ploetzke
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

#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/simulation.h"
#include "ide_secir/parameters_io.h"
#include "memilio/config.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/date.h"
#include <iostream>

int main()
{
    ScalarType N      = 10000;
    ScalarType deaths = 13.10462213;
    ScalarType dt     = 0.5;

    // Initialize model.
    mio::isecir::Model model(mio::TimeSeries<ScalarType>((int)mio::isecir::InfectionTransition::Count), N, deaths);
    // Attention: This example is only working if the file cases_all_germany_ma7.json is previously downloaded and stored in the right folder.
    // TODO: write directions how to dowload and what format.
    auto status = mio::isecir::set_initial_flows(model, dt, "../../data/pydata/Germany/cases_all_germany_ma7.json",
                                                 mio::Date(2020, 12, 24));
    if (!status) {
        std::cout << "Error: " << status.error().formatted_message();
    }
    model.m_transitions.print_table({"S->E", "E->C", "C->I", "C->R", "I->H", "I->R", "H->U", "H->R", "U->D", "U->R"},
                                    16, 8);
    // Carry out simulation.
    /*mio::isecir::Simulation sim(model, 0, dt);

    sim.get_transitions().print_table({"S->E", "E->C", "C->I", "C->R", "I->H", "I->R", "H->U", "H->R", "U->D", "U->R"},
                                      16, 8);*/
}
