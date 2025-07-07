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

int main()
{
    // using Vec = mio::TimeSeries<ScalarType>::Vector;

    // ScalarType tmax = 10.;
    // mio::unused(tmax);
    // ScalarType dt                  = 1;
    // size_t gregory_order           = 1;
    // size_t finite_difference_order = 1; // Possibilities are 1 or 4

    // Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));
    // vec_init[(size_t)mio::isir::InfectionState::Susceptible] = 9910.;
    // // Scheme currently only works if Infected=0 in the beginning.
    // vec_init[(size_t)mio::isir::InfectionState::Infected]  = 0.;
    // vec_init[(size_t)mio::isir::InfectionState::Recovered] = 90.;

    // ScalarType N = vec_init.sum();

    // mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);

    // // TODO: it would be sufficient to have finite_difference_order of time steps before gregory_order
    // init_populations.add_time_point(-(ScalarType)finite_difference_order * dt, vec_init);
    // while (init_populations.get_last_time() < (gregory_order - 1) * dt - 1e-10) {
    //     init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
    // }

    // // Initialize model.
    // mio::isir::Model model(std::move(init_populations), N, gregory_order, finite_difference_order);

    // mio::ExponentialSurvivalFunction exp(1 / 2.);
    // mio::StateAgeFunctionWrapper dist(exp);
    // std::vector<mio::StateAgeFunctionWrapper> vec_dist((size_t)mio::isir::InfectionTransition::Count, dist);
    // model.parameters.get<mio::isir::TransitionDistributions>() = vec_dist;

    // mio::ConstantFunction constfunc(0.1);
    // mio::StateAgeFunctionWrapper const_wrapper(constfunc);
    // model.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = const_wrapper;
    // model.parameters.get<mio::isir::RiskOfInfectionFromSymptomatic>()   = const_wrapper;

    // // Carry out simulation.
    // mio::isir::Simulation sim(model, dt);
    // sim.advance(tmax);

    // sim.get_result().print_table({"S", "I", "R"}, 16, 18);

    // // sim.get_flows().print_table();

    // sim.get_susceptibles_difference().print_table({}, 16, 18);
}
