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

ScalarType t0   = 0.;
ScalarType tmax = 1.;
} // namespace params

mio::IOResult<void> simulate_ide(std::vector<ScalarType> ide_exponents, size_t gregory_order, std::string save_dir = "")
{
    using namespace params;
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    // Only set S because this is the only compartment we consider at the moment.
    Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));
    vec_init[(size_t)mio::isir::InfectionState::Susceptible] = 90.;

    ScalarType total_population = 100.;

    for (ScalarType ide_exponent : ide_exponents) {

        ScalarType dt = pow(10, -ide_exponent);
        std::cout << "Simulate with dt = " << dt << std::endl;

        // Add time points S_0, S_1, S_{n0-1} to init_populations as these values are assumed to be known.
        mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);
        init_populations.add_time_point(0, vec_init);
        while (init_populations.get_last_time() < (gregory_order - 1) * dt - 1e-10) {
            init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
        }

        // init_populations.print_table();

        // Initialize model.
        mio::isir::ModelMessina model(std::move(init_populations), total_population, gregory_order);

        mio::NormalDistributionDensity normaldensity(0.4, 0.6);
        mio::StateAgeFunctionWrapper dist(normaldensity);
        std::vector<mio::StateAgeFunctionWrapper> vec_dist((size_t)mio::isir::InfectionTransition::Count, dist);
        model.parameters.get<mio::isir::TransitionDistributions>() = vec_dist;

        mio::ConstantFunction transmissiononcontact(1.5);
        mio::StateAgeFunctionWrapper transmissiononcontact_wrapper(transmissiononcontact);
        model.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = transmissiononcontact_wrapper;

        mio::ConstantFunction riskofinfection(1.);
        mio::StateAgeFunctionWrapper riskofinfection_wrapper(riskofinfection);
        model.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = riskofinfection_wrapper;

        mio::ContactMatrixGroup contact_matrix             = mio::ContactMatrixGroup(1, 1);
        contact_matrix[0]                                  = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 0.001));
        model.parameters.get<mio::isir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

        // Carry out simulation.
        mio::isir::SimulationMessina sim(model, dt);
        sim.advance_messina(tmax);

        // sim.get_result().print_table();

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

    return mio::success();
}

int main()
{
    std::string result_dir = "../../simulation_results/messina_new/";
    // Make folder if not existent yet.
    boost::filesystem::path dir(result_dir);
    boost::filesystem::create_directories(dir);

    // std::vector<size_t> gregory_orders    = {1, 2, 3};
    // std::vector<ScalarType> ide_exponents = {1, 2, 3};

    // for (size_t gregory_order : gregory_orders) {

    //     std::cout << "Using Gregory order = " << gregory_order << std::endl;
    //     mio::IOResult<void> result = simulate_ide(ide_exponents, gregory_order, result_dir);
    // }

    // // Compute groundtruth.

    size_t gregory_order                              = 3;
    std::vector<ScalarType> ide_exponents_groundtruth = {5};

    std::cout << "Using Gregory order = " << gregory_order << std::endl;
    mio::IOResult<void> result = simulate_ide(ide_exponents_groundtruth, gregory_order, result_dir);
}
