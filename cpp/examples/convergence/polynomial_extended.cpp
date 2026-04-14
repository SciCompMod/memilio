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

#include "ide_sir/model_polynomial_extended.h"
#include "ide_sir/infection_state.h"
#include "ide_sir/simulation_polynomial_extended.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/io/result_io.h"
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <vector>

using namespace mio;
namespace params
{
size_t num_agegroups = 1;

ScalarType TransmissionProbabilityOnContact = 1.;
ScalarType RiskOfInfectionFromSymptomatic   = 1.;
ScalarType Seasonality                      = 0.;

ScalarType cont_freq = 0.7;

ScalarType S0 = 999000.;
ScalarType I0 = 1000.;
ScalarType R0 = 0.;
// ScalarType total_population = S0 + I0 + R0;
ScalarType total_population = 100;

// ScalarType t_init = 0.;
} // namespace params

ScalarType polynomial(ScalarType time, ScalarType alpha, ScalarType t_init)
{
    using namespace params;

    return alpha / 30. * std::pow(time, 6) + (1 - alpha * time) / 5. * std::pow(time, 5) -
           (1 - alpha * (time - t_init)) * total_population * time - alpha * total_population / 2. * std::pow(time, 2);
}

mio::IOResult<void> simulate_ide(std::vector<ScalarType> ide_exponents, size_t gregory_order, ScalarType t0,
                                 ScalarType t_init, ScalarType tmax, ScalarType alpha, bool single_integral,
                                 std::string save_dir = "")
{
    using namespace params;
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    for (ScalarType ide_exponent : ide_exponents) {

        ScalarType dt_ide = pow(10, -ide_exponent);
        std::cout << "Simulation with dt=" << dt_ide << std::endl;

        std::cout << "Initializing with given groundtruth.\n";

        // Initialize time points before t0_ide based on groundtruth.

        Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));

        // Add values to init_populations.
        vec_init[(size_t)mio::isir::InfectionState::Susceptible] = polynomial(t0, alpha, t_init);

        mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);
        init_populations.add_time_point(t0, vec_init);

        while (init_populations.get_last_time() < t_init - 1e-10) {
            ScalarType time = init_populations.get_last_time() + dt_ide;
            // vec_init[(size_t)mio::isir::InfectionState::Susceptible] = std::pow(time, 3) + std::pow(time, 2);
            // vec_init[(size_t)mio::isir::InfectionState::Susceptible] = std::pow(time, 3) + std::pow(time, 2);
            vec_init[(size_t)mio::isir::InfectionState::Susceptible] = polynomial(time, alpha, t_init);
            init_populations.add_time_point(time, vec_init);
        }

        mio::TimeSeries<ScalarType> init_populations_copy = init_populations;

        // Initialize model.
        mio::isir::ModelPolynomialExtended model(std::move(init_populations), gregory_order, total_population);

        // Carry out simulation.
        mio::isir::SimulationPolynomialExtended sim(model, dt_ide);
        if (single_integral) {
            sim.advance_single(tmax, alpha);
        }
        else {
            sim.advance_double(tmax, alpha);
        }

        if (!save_dir.empty()) {
            // Save compartments.
            mio::TimeSeries<ScalarType> compartments = sim.get_result();
            auto save_result_status_ide =
                mio::save_result({compartments}, {0}, num_agegroups,
                                 save_dir + "result_ide_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
                                     "_gregoryorder=" + fmt::format("{}", gregory_order) + ".h5");
        }
    }

    return mio::success();
}

int main()
{
    /* In this example we want to examine the convergence behavior under the assumption of exponential stay time 
    distributions. In this case, we can compare the solution of the IDE simulation with a corresponding ODE solution. */

    using namespace params;

    ScalarType t0     = 0.;
    ScalarType t_init = 2.;
    ScalarType tmax   = 4.;

    ScalarType alpha = 0.5;

    bool single_integral = false;

    std::vector<ScalarType> ide_exponents = {0, 1, 2};
    std::vector<size_t> gregory_orders    = {1, 2, 3};

    std::string save_dir =
        fmt::format("../../simulation_results/2026-04-14/convergence_polynomial_ext_alpha={}_singleint={}/"
                    "polynomial_t0={}_tinit={}_tmax={}/",
                    alpha, single_integral, t0, t_init, tmax);

    // Make folder if not existent yet.
    boost::filesystem::path dir(save_dir);
    boost::filesystem::create_directories(dir);

    // Do IDE simulations.
    for (size_t gregory_order : gregory_orders) {
        std::cout << std::endl;
        std::cout << "Gregory order: " << gregory_order << std::endl;
        mio::IOResult<void> result_ide =
            simulate_ide(ide_exponents, gregory_order, t0, t_init, tmax, alpha, single_integral, save_dir);
    }
}
// }
