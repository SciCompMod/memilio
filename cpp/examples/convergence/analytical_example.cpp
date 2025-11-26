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

#include "ide_sir_analytical/model.h"
#include "ide_sir_analytical/simulation.h"
#include "ide_sir_analytical/infection_state.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/io/result_io.h"
#include <vector>
#include <math.h>

using namespace mio;
namespace params
{

size_t num_agegroups = 1;

ScalarType y0 = 1.;

} // namespace params

mio::IOResult<void> simulate_ide(std::vector<ScalarType> ide_exponents, size_t gregory_order, ScalarType t0_ide,
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
        vec_init[(size_t)mio::isir::InfectionState::Susceptible] = cosh(0.);
        vec_init[(size_t)mio::isir::InfectionState::Infected]    = 0.;
        vec_init[(size_t)mio::isir::InfectionState::Recovered]   = 0.;

        init_populations.add_time_point(t0_ide, vec_init);

        while (init_populations.get_last_time() < 3 * dt_ide - 1e-10) {
            init_populations.add_time_point(init_populations.get_last_time() + dt_ide);
            // vec_init[(size_t)mio::isir::InfectionState::Susceptible] =
            //     cosh(pow(init_populations.get_last_time(), 2) / 2.);
            vec_init[(size_t)mio::isir::InfectionState::Susceptible] = cosh(init_populations.get_last_time());
            vec_init[(size_t)mio::isir::InfectionState::Infected]    = 0.;
            vec_init[(size_t)mio::isir::InfectionState::Recovered]   = 0.;
            init_populations.get_last_value()                        = vec_init;
        }

        // Initialize model.
        mio::isir::ModelAnalytical model(std::move(init_populations), gregory_order);

        // Carry out simulation.
        mio::isir::SimulationAnalytical sim(model, dt_ide);
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
    /* In this example we want to examine the convergence behavior under the assumption of exponential stay time 
    distributions. In this case, we can compare the solution of the IDE simulation with a corresponding ODE solution. */

    using namespace params;

    ScalarType t0   = 0.;
    ScalarType tmax = 5.;

    std::vector<ScalarType> ide_exponents = {0, 1, 2, 3, 4};
    std::vector<size_t> gregory_orders    = {1, 2, 3};

    std::string save_dir = fmt::format("../../simulation_results/2025-11-17/analytical_example/"
                                       "t0ide={}_tmax={}/",
                                       t0, tmax);

    // Make folder if not existent yet.
    boost::filesystem::path dir(save_dir);
    boost::filesystem::create_directories(dir);

    // Do IDE simulations.
    for (size_t gregory_order : gregory_orders) {
        std::cout << std::endl;
        std::cout << "Gregory order: " << gregory_order << std::endl;
        mio::IOResult<void> result_ide = simulate_ide(ide_exponents, gregory_order, t0, tmax, save_dir);
    }
}