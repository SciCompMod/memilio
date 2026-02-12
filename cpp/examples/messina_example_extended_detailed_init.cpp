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
#include "ode_sir/model.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/utils/time_series.h"
#include "memilio/io/result_io.h"
#include <vector>

namespace params
{
size_t num_agegroups = 1;

// Start time of initial values
ScalarType t_init = -0.5;
// Simulation start
ScalarType t0 = 0.;
// Simulation end
ScalarType tmax = 1.;

ScalarType S0 = 90.;
ScalarType I0 = 10.;
ScalarType R0 = 0.;

ScalarType total_population = S0 + I0 + R0;
ScalarType beta             = 0.001;
ScalarType cont_freq        = beta * total_population;

size_t finite_difference_order = 4;
} // namespace params

// This function returns a TimeSeries containing the groundtruth, which is either the one given to this function or the
// one computed here.
mio::IOResult<mio::TimeSeries<ScalarType>>
simulate_ide(std::vector<ScalarType> ide_exponents, size_t gregory_order, std::string save_dir = "",
             mio::TimeSeries<ScalarType> result_groundtruth =
                 mio::TimeSeries<ScalarType>((size_t)mio::isir::InfectionState::Count))
{
    using namespace params;
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    for (ScalarType ide_exponent : ide_exponents) {
        ScalarType dt = pow(10., -ide_exponent);
        std::cout << "Simulate with dt = " << dt << std::endl;

        // Define init_populations; depending on whether we have a groundtruth available or if we are computing it here,
        // we are setting init_populations differently.
        mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);

        // std::cout << "Num time points groundtruth: " << result_groundtruth.get_num_time_points() << std::endl;

        // Check if groundtruth is available and initialize accordingly.
        if (result_groundtruth.get_num_time_points() == 0) {
            // Initialize first (-t_init/dt)+(gregory_order-1) time points with constant values.
            Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));
            vec_init[(size_t)mio::isir::InfectionState::Susceptible] = S0;
            vec_init[(size_t)mio::isir::InfectionState::Infected]    = I0;
            vec_init[(size_t)mio::isir::InfectionState::Recovered]   = R0;
            // Add time points from t_init to t0 to init_populations.
            init_populations.add_time_point(t_init, vec_init);
            while (init_populations.get_last_time() < t0 * dt - 1e-10) {
                init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
            }
        }
        else {
            // Initialize first (gregory_order-1) time points based on groundtruth.
            ScalarType dt_groundtruth = result_groundtruth.get_time(1) - result_groundtruth.get_time(0);

            size_t groundtruth_index = size_t(dt / dt_groundtruth);

            Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));

            std::vector<size_t> compartments = {(size_t)mio::isir::InfectionState::Susceptible,
                                                (size_t)mio::isir::InfectionState::Infected,
                                                (size_t)mio::isir::InfectionState::Recovered};

            // Add values for t=0.
            for (size_t compartment : compartments) {
                vec_init[compartment] = result_groundtruth.get_value(0)[compartment];
            }
            init_populations.add_time_point(t_init, vec_init);

            // Add values for t_init,...,t0..
            while (init_populations.get_last_time() < t0 - 1e-10) {
                for (size_t compartment : compartments) {
                    vec_init[compartment] = result_groundtruth.get_value(
                        size_t(init_populations.get_num_time_points() * groundtruth_index))[compartment];
                }
                init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
            }
        }

        // Initialize model.
        mio::isir::ModelMessinaExtendedDetailedInit model(std::move(init_populations), total_population, gregory_order,
                                                          finite_difference_order);

        mio::NormalDistributionDensity normaldensity(0.4, 0.6);
        mio::StateAgeFunctionWrapper dist(normaldensity);
        // mio::ExponentialSurvivalFunction exp(2.);
        // mio::StateAgeFunctionWrapper dist(exp);
        std::vector<mio::StateAgeFunctionWrapper<ScalarType>> vec_dist((size_t)mio::isir::InfectionTransition::Count,
                                                                       dist);
        model.parameters.get<mio::isir::TransitionDistributions>() = vec_dist;

        mio::ConstantFunction transmissiononcontact(1.5);
        mio::StateAgeFunctionWrapper transmissiononcontact_wrapper(transmissiononcontact);
        model.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = transmissiononcontact_wrapper;

        mio::ConstantFunction riskofinfection(1.);
        mio::StateAgeFunctionWrapper riskofinfection_wrapper(riskofinfection);
        model.parameters.get<mio::isir::RiskOfInfectionFromSymptomatic>() = riskofinfection_wrapper;

        mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup<ScalarType>(1, 1);
        contact_matrix[0] =
            mio::ContactMatrix<ScalarType>(Eigen::MatrixXd::Constant(1, 1, cont_freq / total_population));
        model.parameters.get<mio::isir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

        // Carry out simulation.
        mio::isir::SimulationMessinaExtendedDetailedInit sim(model, dt);
        sim.advance(tmax);

        // If no groundtruth is gibven as input, we set the here computed results as groundtruth.
        if (result_groundtruth.get_num_time_points() == 0) {
            result_groundtruth = sim.get_result();
        }

        // sim.get_result().print_table();

        mio::unused(save_dir);

        if (!save_dir.empty()) {
            // Save compartments.
            mio::TimeSeries<ScalarType> compartments = sim.get_result();
            auto save_result_status_ide =
                mio::save_result({compartments}, {0}, num_agegroups,
                                 save_dir + "result_ide_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
                                     "_gregoryorder=" + fmt::format("{}", gregory_order) + ".h5");

            if (!save_result_status_ide) {
                return mio::failure(mio::StatusCode::InvalidValue,
                                    "Error occured while saving the IDE simulation results.");
            }
        }
    }

    return mio::success(result_groundtruth);
}

int main()
{
    std::string result_dir = "../../simulation_results/messina_model_extended_detailed_init_test/";
    // Make folder if not existent yet.
    boost::filesystem::path dir(result_dir);
    boost::filesystem::create_directories(dir);

    // Compute groundtruth.

    size_t gregory_order_groundtruth                  = 3;
    std::vector<ScalarType> ide_exponents_groundtruth = {4};

    std::cout << "Using Gregory order = " << gregory_order_groundtruth << std::endl;
    mio::IOResult<mio::TimeSeries<ScalarType>> result_groundtruth =
        simulate_ide(ide_exponents_groundtruth, gregory_order_groundtruth, result_dir);

    // // Simulate with larger timesteps and use result_groundtruth for initialization.
    // std::vector<size_t> gregory_orders    = {1, 2, 3};
    // std::vector<ScalarType> ide_exponents = {1, 2, 3, 4};

    // for (size_t gregory_order : gregory_orders) {
    //     std::cout << "Using Gregory order = " << gregory_order << std::endl;
    //     mio::IOResult<mio::TimeSeries<ScalarType>> result =
    //         simulate_ide(ide_exponents, gregory_order, result_dir, result_groundtruth.value());
    // }
}
