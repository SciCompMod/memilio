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
ScalarType tmax = 10.;

ScalarType TimeInfected                     = 1.4;
ScalarType TransmissionProbabilityOnContact = 0.1;
ScalarType RiskOfInfectionFromSymptomatic   = 0.1;
ScalarType Seasonality                      = 0.;
ScalarType cont_freq                        = 1.;
} // namespace params

// void simulate_ode(ScalarType dt)
// {
//     using namespace params;

//     mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

//     mio::osir::Model<ScalarType> model(num_agegroups);

//     ScalarType total_population = 10000;

//     model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Infected}]  = 1000;
//     model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Recovered}] = 1000;
//     model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Susceptible}] =
//         total_population - model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Infected}] -
//         model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Recovered}];
//     model.parameters.set<mio::osir::TimeInfected<ScalarType>>(2);
//     model.parameters.set<mio::osir::TransmissionProbabilityOnContact<ScalarType>>(0.5);

//     mio::ContactMatrixGroup& contact_matrix =
//         model.parameters.get<mio::osir::ContactPatterns<ScalarType>>().get_cont_freq_mat();
//     contact_matrix[0].get_baseline().setConstant(2.7);
//     contact_matrix[0].add_damping(0.6, mio::SimulationTime(12.5));
//     model.check_constraints();

//     std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator =
//         std::make_shared<mio::EulerIntegratorCore<ScalarType>>();
//     auto sir = simulate(t0, tmax, dt, model, integrator);
// }

mio::IOResult<void> simulate_ide(std::vector<ScalarType> ide_exponents, size_t gregory_order,
                                 size_t finite_difference_order, std::string save_dir = "")
{
    using namespace params;
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));
    vec_init[(size_t)mio::isir::InfectionState::Susceptible] = 9910.;
    // Scheme currently only works if Infected=0 in the beginning.
    vec_init[(size_t)mio::isir::InfectionState::Infected]  = 0.;
    vec_init[(size_t)mio::isir::InfectionState::Recovered] = 90.;

    ScalarType total_population = vec_init.sum();

    for (ScalarType ide_exponent : ide_exponents) {

        ScalarType dt = pow(10, -ide_exponent);
        std::cout << "Simulation with " << dt << std::endl;

        // TODO: it would be sufficient to have finite_difference_order of time steps before gregory_order
        mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);
        init_populations.add_time_point(-(ScalarType)finite_difference_order * dt, vec_init);
        while (init_populations.get_last_time() < (gregory_order - 1) * dt - 1e-10) {
            init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
        }

        // Initialize model.
        mio::isir::Model model(std::move(init_populations), total_population, gregory_order, finite_difference_order);

        mio::ExponentialSurvivalFunction exp(1. / TimeInfected);
        mio::StateAgeFunctionWrapper dist(exp);
        std::vector<mio::StateAgeFunctionWrapper> vec_dist((size_t)mio::isir::InfectionTransition::Count, dist);
        model.parameters.get<mio::isir::TransitionDistributions>() = vec_dist;

        mio::ConstantFunction transmissiononcontact(TransmissionProbabilityOnContact);
        mio::StateAgeFunctionWrapper transmissiononcontact_wrapper(transmissiononcontact);
        model.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = transmissiononcontact_wrapper;

        mio::ConstantFunction riskofinfection(RiskOfInfectionFromSymptomatic);
        mio::StateAgeFunctionWrapper riskofinfection_wrapper(transmissiononcontact);
        model.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = riskofinfection_wrapper;

        // Carry out simulation.
        mio::isir::Simulation sim(model, dt);
        sim.advance(tmax);

        if (!save_dir.empty()) {
            // Save compartments.
            mio::TimeSeries<ScalarType> compartments = sim.get_result();
            auto save_result_status_ide =
                mio::save_result({compartments}, {0}, num_agegroups,
                                 save_dir + "result_ide_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
                                     "_gregoryorder=" + fmt::format("{}", gregory_order) +
                                     "_finitedifforder=" + fmt::format("{}", finite_difference_order) + ".h5");

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
    std::string result_dir = "../../simulation_results/";
    // Make folder if not existent yet.
    boost::filesystem::path dir(result_dir);
    boost::filesystem::create_directories(dir);

    size_t finite_difference_order = 1;

    std::vector<ScalarType> ide_exponents = {0, 1, 2};

    size_t gregory_order = 1;

    mio::IOResult<void> result = simulate_ide(ide_exponents, gregory_order, finite_difference_order, result_dir);

    // Get groundtruth with gregory_order = 3
    gregory_order = 3;
    ide_exponents = {3};

    result = simulate_ide(ide_exponents, gregory_order, finite_difference_order, result_dir);
}
