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
#include "ide_secir/initialflows.h"
#include "memilio/config.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"

int main()
{
    ScalarType N      = 10000;
    ScalarType deaths = 13.10462213;
    ScalarType dt     = 1;

    // Initialize model.
    mio::isecir::Model model(mio::TimeSeries<ScalarType>((int)mio::isecir::InfectionTransition::Count), N, deaths);

    ScalarType rki_cases_dummy{10.};
    ScalarType rki_deaths_dummy{2.};
    mio::isecir::set_initial_flows(model, dt, rki_cases_dummy, rki_deaths_dummy);

    // Carry out simulation.
    mio::isecir::Simulation sim(model, 0, dt);

    sim.get_transitions().print_table({"S->E", "E->C", "C->I", "C->R", "I->H", "I->R", "H->U", "H->R", "U->D", "U->R"},
                                      16, 8);
}
