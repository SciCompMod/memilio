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
#include "ide_sir/model_polynomial.h"
#include "ide_sir/infection_state.h"
#include "ide_sir/simulation_polynomial.h"
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

ScalarType S0               = 999000.;
ScalarType I0               = 1000.;
ScalarType R0               = 0.;
ScalarType total_population = S0 + I0 + R0;
} // namespace params

mio::IOResult<void> simulate_ide(std::vector<ScalarType> ide_exponents, size_t gregory_order, ScalarType t0,
                                 ScalarType t_init, ScalarType tmax, std::string save_dir = "")
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
        vec_init[(size_t)mio::isir::InfectionState::Susceptible] = 0.;

        mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);
        init_populations.add_time_point(t0, vec_init);

        while (init_populations.get_last_time() < t_init - 1e-10) {
            ScalarType time = init_populations.get_last_time() + dt_ide;
            // vec_init[(size_t)mio::isir::InfectionState::Susceptible] = std::pow(time, 3) + std::pow(time, 2);
            // vec_init[(size_t)mio::isir::InfectionState::Susceptible] = std::pow(time, 3) + std::pow(time, 2);
            vec_init[(size_t)mio::isir::InfectionState::Susceptible] = std::pow(time, 5);
            init_populations.add_time_point(time, vec_init);
        }

        mio::TimeSeries<ScalarType> init_populations_copy = init_populations;

        // Initialize model.
        mio::isir::ModelPolynomial model_single(std::move(init_populations), gregory_order);
        mio::isir::ModelPolynomial model_double(std::move(init_populations_copy), gregory_order);

        // Carry out simulation.
        // mio::isir::SimulationPolynomial sim_single(model_single, dt_ide);
        // sim_single.advance_single(tmax);

        mio::isir::SimulationPolynomial sim_double(model_double, dt_ide);
        sim_double.advance_double(tmax);

        if (!save_dir.empty()) {
            // Save compartments.
            // mio::TimeSeries<ScalarType> compartments_single = sim_single.get_result();
            // auto save_result_status_ide =
            //     mio::save_result({compartments_single}, {0}, num_agegroups,
            //                      save_dir + "result_ide_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
            //                          "_gregoryorder=" + fmt::format("{}", gregory_order) + ".h5");

            mio::TimeSeries<ScalarType> compartments_double = sim_double.get_result();
            auto save_result_status_ide_double =
                mio::save_result({compartments_double}, {0}, num_agegroups,
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

    ScalarType t0                         = 0.;
    std::vector<ScalarType> t_init_values = {5.};
    ScalarType tmax                       = 10.;

    std::vector<ScalarType> ide_exponents = {0, 1, 2};
    std::vector<size_t> gregory_orders    = {1, 2, 3};

    // bool single_integral = true;

    for (ScalarType t_init : t_init_values) {

        // for (size_t num_days : num_days_vec) {
        // ScalarType tmax = t0_ide + num_days;

        std::string save_dir =
            fmt::format("../../simulation_results/2026-04-09/convergence_polynomial_order5/double_integral/"
                        "polynomial_t0={}_tinit={}_tmax={}/",
                        t0, t_init, tmax);

        // Make folder if not existent yet.
        boost::filesystem::path dir(save_dir);
        boost::filesystem::create_directories(dir);

        // Do IDE simulations.
        for (size_t gregory_order : gregory_orders) {
            std::cout << std::endl;
            std::cout << "Gregory order: " << gregory_order << std::endl;
            mio::IOResult<void> result_ide = simulate_ide(ide_exponents, gregory_order, t0, t_init, tmax, save_dir);
        }
    }
}
