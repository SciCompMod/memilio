/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler
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

#include "ide_sir/model.h"
#include "ide_sir/infection_state.h"
#include "ide_sir/parameters.h"
#include "ide_sir/simulation.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/utils/time_series.h"
#include "ode_secir/parameters.h"

int main()
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax      = 4.;
    ScalarType dt        = 1.;
    size_t gregory_order = 2;

    Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));
    vec_init[(size_t)mio::isir::InfectionState::Susceptible] = 900000;
    vec_init[(size_t)mio::isir::InfectionState::Infected]    = 1;
    vec_init[(size_t)mio::isir::InfectionState::Recovered]   = 900;

    ScalarType N = vec_init.sum();

    mio::TimeSeries<ScalarType> init_populations((int)mio::isir::InfectionState::Count);

    init_populations.add_time_point(0., vec_init);
    while (init_populations.get_last_time() < (gregory_order - 1) * dt - 1e-6) {
        init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isir::Model model(std::move(init_populations), N);

    mio::ExponentialSurvivalFunction exp(9.);
    mio::StateAgeFunctionWrapper dist(exp);
    std::vector<mio::StateAgeFunctionWrapper> vec_dist((size_t)mio::isir::InfectionTransition::Count, dist);
    model.parameters.get<mio::isir::TransitionDistributions>() = vec_dist;

    mio::ConstantFunction constfunc(0.1);
    mio::StateAgeFunctionWrapper const_wrapper(constfunc);
    model.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = const_wrapper;

    // Carry out simulation.
    mio::isir::Simulation sim(model, dt, gregory_order);
    sim.advance(tmax);

    sim.get_result().print_table({"S", "I", "R"}, 16, 8);
}
