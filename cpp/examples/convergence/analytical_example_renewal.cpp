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

#include "ide_sir/infection_state.h"
#include "ide_sir/parameters.h"
#include "ide_sir_analytical_renewal/model.h"
#include "ide_sir_analytical_renewal/simulation.h"
#include "ide_sir_analytical_renewal/infection_state.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/utils/time_series.h"
#include "memilio/io/result_io.h"
#include <vector>
#include <math.h>

using namespace mio;
namespace params
{
size_t num_agegroups = 1;

size_t finite_difference_order = 4;

ScalarType t0_init = 10.;

ScalarType rho   = 5.;
ScalarType kappa = 5.1;
} // namespace params

mio::IOResult<void> simulate_ide(std::vector<ScalarType> ide_exponents, size_t gregory_order, ScalarType t0_groundtruth,
                                 ScalarType tmax, std::string save_dir = "")
{
    using namespace params;
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    for (ScalarType ide_exponent : ide_exponents) {

        ScalarType dt_ide = pow(10, -ide_exponent);
        std::cout << "Simulation with dt=" << dt_ide << std::endl;

        mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);

        std::cout << "Initializing with given groundtruth.\n";

        Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));

        // Add values to init_populations.
        vec_init[(size_t)mio::isir::InfectionState::Susceptible] = rho;
        vec_init[(size_t)mio::isir::InfectionState::Infected]    = 0.;
        vec_init[(size_t)mio::isir::InfectionState::Recovered]   = 0.;

        init_populations.add_time_point(t0_groundtruth, vec_init);

        while (init_populations.get_last_time() < t0_init - 1e-10) {
            ScalarType time = init_populations.get_last_time() + dt_ide;

            vec_init[(size_t)mio::isir::InfectionState::Susceptible] = rho * std::exp((rho - kappa) * time);
            vec_init[(size_t)mio::isir::InfectionState::Infected] =
                (rho - kappa) * (std::exp((rho - kappa) * time) - std::exp(-kappa * time));
            vec_init[(size_t)mio::isir::InfectionState::Recovered] =
                vec_init[(size_t)mio::isir::InfectionState::Susceptible] - rho -
                vec_init[(size_t)mio::isir::InfectionState::Infected];
            init_populations.add_time_point(time, vec_init);
        }

        // Initialize model.
        mio::isir::ModelAnalyticalRenewal model(std::move(init_populations), gregory_order, finite_difference_order);

        mio::ConstantFunction<ScalarType> rho_func(rho);
        mio::StateAgeFunctionWrapper<ScalarType> saf_wrapper(rho_func);
        model.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = saf_wrapper;

        mio::ExponentialSurvivalFunction exp(kappa);
        mio::StateAgeFunctionWrapper dist(exp);
        std::vector<mio::StateAgeFunctionWrapper<ScalarType>> vec_dist((size_t)mio::isir::InfectionTransition::Count,
                                                                       dist);
        model.parameters.get<mio::isir::TransitionDistributions>() = vec_dist;

        // Carry out simulation.
        mio::isir::SimulationAnalyticalRenewal sim(model, dt_ide);
        sim.advance(tmax);

        if (!save_dir.empty()) {
            // Save compartments.
            mio::TimeSeries<ScalarType> result = sim.get_result();
            auto save_result_status_ide =
                mio::save_result({result}, {0}, num_agegroups,
                                 save_dir + "result_ide_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
                                     "_gregoryorder=" + fmt::format("{}", gregory_order) + ".h5");

            if (!save_result_status_ide) {
                return mio::failure(mio::StatusCode::InvalidValue,
                                    "Error occured while saving the IDE simulation results.");
            }
        }
    }

    return mio::success();
}

int main()
{
    using namespace params;

    ScalarType t0   = 0.;
    ScalarType tmax = 50.;

    std::vector<ScalarType> ide_exponents = {0, 1};
    std::vector<size_t> gregory_orders    = {3};

    std::string save_dir = fmt::format("../../simulation_results/2026-06-14/analytical_renewal_rho={}_kappa={}/"
                                       "t0ide={}_tinit={}_tmax={}/",
                                       rho, kappa, t0, t0_init, tmax);

    std::cout << save_dir << std::endl;

    // Make folder if not existent yet.
    std::filesystem::path dir(save_dir);
    std::filesystem::create_directories(dir);

    // Do IDE simulations.
    for (size_t gregory_order : gregory_orders) {
        std::cout << std::endl;
        std::cout << "Gregory order: " << gregory_order << std::endl;
        mio::IOResult<void> result_ide = simulate_ide(ide_exponents, gregory_order, t0, tmax, save_dir);
    }
}