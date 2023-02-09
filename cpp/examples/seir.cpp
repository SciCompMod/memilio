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
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "memilio/math/integrator.h"
#include "memilio/math/stepper_wrapper.h"
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/math/flow_calculator.h"
#include <iostream>
#include <type_traits>

void print_to_terminal(const mio::TimeSeries<ScalarType>& results, const std::vector<std::string>& state_names)
{
    // print column labels
    printf("%-16s  ", "Time");
    for (size_t k = 0; k < static_cast<size_t>(results.get_num_elements()); k++) {
        if (k < state_names.size()) {
            printf(" %-16s", state_names[k].data()); // print underlying char*
        }
        else {
            printf(" %-16s", ("#" + std::to_string(k + 1)).data());
        }
    }
    // print values as table
    auto num_points = static_cast<size_t>(results.get_num_time_points());
    for (size_t i = 0; i < num_points; i++) {
        printf("\n%16.6f", results.get_time(i));
        auto res_i = results.get_value(i);
        for (size_t j = 0; j < static_cast<size_t>(res_i.size()); j++) {
            printf(" %16.6f", res_i[j]);
        }
    }
    printf("\n");
}

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 10;
    double dt   = 0.001;

    mio::log_info("Simulating SEIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::oseir::Model model;

    double total_population                                                                            = 10000;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}]   = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}]  = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}] = 100;
    model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Exposed)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::oseir::InfectionState>(mio::oseir::InfectionState::Recovered)}];
    // suscetible now set with every other update
    // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
    model.parameters.set<mio::oseir::TimeExposed>(5.2);
    model.parameters.set<mio::oseir::TimeInfected>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact>(0.04);
    model.parameters.get<mio::oseir::ContactPatterns>().get_baseline()(0, 0) = 10;

    model.check_constraints();
    // print_seir_params(model);
    /* auto I = std::make_shared<mio::flow_calculator>(
        model.parameters.get<mio::oseir::Flows>(),
        mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>()); */
    auto seir = simulate(t0, tmax, dt, model);

    print_to_terminal(seir, {});

    printf("\n number total: %f\n",
           seir.get_last_value()[0] + seir.get_last_value()[1] + seir.get_last_value()[2] + seir.get_last_value()[3]);
}
