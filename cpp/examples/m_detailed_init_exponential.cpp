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
#include "memilio/utils/logging.h"
#include "ode_sir/model.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/utils/time_series.h"
#include "memilio/io/result_io.h"
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <vector>

namespace params
{
size_t num_agegroups = 1;

ScalarType t0_ode = 0.;
ScalarType t0_ide = 50.;
ScalarType tmax   = 51.;

ScalarType TimeInfected = 2.;
// This parameter is chosen differently than in the example from the paper, as this is not a valid choice for a probability.
// Instead we scale the contact frequency with a factor of 1.5.
ScalarType TransmissionProbabilityOnContact = 1.;
ScalarType RiskOfInfectionFromSymptomatic   = 1.;
ScalarType Seasonality                      = 0.;

ScalarType S0 = 95.;
ScalarType I0 = 5.;
ScalarType R0 = 0.;

ScalarType total_population = S0 + I0 + R0;
// Note that the contacts are currently differently defined in ODE and IDE model, which is why they are set differently
// according to the contact frequency.
ScalarType cont_freq = 1.5 * 0.1;
ScalarType beta      = cont_freq / total_population;
} // namespace params

mio::IOResult<mio::TimeSeries<ScalarType>> simulate_ode(ScalarType ode_exponent, std::string save_dir = "")
{
    using namespace params;

    ScalarType dt_ode = pow(10, -ode_exponent);

    mio::log_info("Simulating ODE-SIR; t={} ... {} with dt = {}.", t0_ode, tmax, dt_ode);

    mio::osir::Model<ScalarType> model(num_agegroups);

    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Susceptible}] = S0;
    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Infected}]    = I0;
    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Recovered}]   = R0;

    model.parameters.set<mio::osir::TimeInfected<ScalarType>>(TimeInfected);
    model.parameters.set<mio::osir::TransmissionProbabilityOnContact<ScalarType>>(TransmissionProbabilityOnContact);

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    model.parameters.get<mio::osir::ContactPatterns<ScalarType>>() = mio::UncertainContactMatrix(contact_matrix);

    model.check_constraints();

    std::shared_ptr<mio::OdeIntegratorCore<ScalarType>> integrator =
        std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_fehlberg78>>();

    auto sir = simulate<ScalarType, mio::osir::Model<ScalarType>>(t0_ode, tmax, dt_ode, model, integrator);

    if (!save_dir.empty()) {
        // Save compartments.
        mio::TimeSeries<ScalarType> compartments = sir;
        auto save_result_status_ode =
            mio::save_result({compartments}, {0}, num_agegroups,
                             save_dir + "result_ode_dt=1e-" + fmt::format("{:.0f}", ode_exponent) + ".h5");

        if (!save_result_status_ode) {
            return mio::failure(mio::StatusCode::InvalidValue,
                                "Error occured while saving the ODE simulation results.");
        }
    }
    return mio::success(sir);
}

mio::IOResult<void> simulate_ide(std::vector<ScalarType> ide_exponents, size_t gregory_order,
                                 size_t finite_difference_order, std::string save_dir = "",
                                 mio::TimeSeries<ScalarType> result_groundtruth =
                                     mio::TimeSeries<ScalarType>((size_t)mio::isir::InfectionState::Count),
                                 bool backwarts_fd = true)
{
    using namespace params;
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    for (ScalarType ide_exponent : ide_exponents) {

        ScalarType dt_ide = pow(10, -ide_exponent);
        std::cout << "Simulation with " << dt_ide << std::endl;

        mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);

        if (result_groundtruth.get_num_time_points() == 0) {
            std::cout << "No groundtruth was given.\n";
        }
        else {
            std::cout << "Initializing with given groundtruth.\n";

            // Initialize time points before t0_ide based on groundtruth.
            ScalarType dt_groundtruth = result_groundtruth.get_time(1) - result_groundtruth.get_time(0);
            // std::cout << "dt groundtruth: " << dt_groundtruth << std::endl;
            size_t groundtruth_index_factor = size_t(dt_ide / dt_groundtruth);
            // std::cout << "groundtruth_index: " << groundtruth_index << std::endl;

            Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));

            std::vector<size_t> compartments = {(size_t)mio::isir::InfectionState::Susceptible,
                                                (size_t)mio::isir::InfectionState::Infected,
                                                (size_t)mio::isir::InfectionState::Recovered};

            // Add values to init_populations.
            for (size_t compartment : compartments) {
                vec_init[compartment] = result_groundtruth.get_value(0)[compartment];
            }
            init_populations.add_time_point(t0_ode, vec_init);

            while (init_populations.get_last_time() < t0_ide - 1e-10) {
                for (size_t compartment : compartments) {
                    vec_init[compartment] = result_groundtruth.get_value(
                        size_t(init_populations.get_num_time_points() * groundtruth_index_factor))[compartment];
                }
                init_populations.add_time_point(init_populations.get_last_time() + dt_ide, vec_init);
            }
        }

        // Initialize model.
        mio::isir::ModelMessinaExtendedDetailedInit model(std::move(init_populations), total_population, gregory_order,
                                                          finite_difference_order);

        mio::ExponentialSurvivalFunction exp(1. / TimeInfected);
        // std::cout << "max support: " << exp.get_support_max(dt_ide) << std::endl;
        ;
        mio::StateAgeFunctionWrapper dist(exp);
        std::vector<mio::StateAgeFunctionWrapper> vec_dist((size_t)mio::isir::InfectionTransition::Count, dist);
        model.parameters.get<mio::isir::TransitionDistributions>() = vec_dist;

        mio::ConstantFunction transmissiononcontact(TransmissionProbabilityOnContact);
        mio::StateAgeFunctionWrapper transmissiononcontact_wrapper(transmissiononcontact);
        model.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = transmissiononcontact_wrapper;

        mio::ConstantFunction riskofinfection(RiskOfInfectionFromSymptomatic);
        mio::StateAgeFunctionWrapper riskofinfection_wrapper(riskofinfection);
        model.parameters.get<mio::isir::RiskOfInfectionFromSymptomatic>() = riskofinfection_wrapper;

        mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
        contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
        model.parameters.get<mio::isir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

        // Carry out simulation.
        mio::isir::SimulationMessinaExtendedDetailedInit sim(model, dt_ide);
        sim.advance(tmax, backwarts_fd);

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
    // Compute groundtruth with ODE model.
    ScalarType ode_exponent = 6;

    // true means that a backwards_fd scheme is used, false means that a central fd scheme is used
    std::vector<bool> fd_schemes = {true};

    for (bool backwards_fd : fd_schemes) {

        std::vector<size_t> finite_difference_orders;
        if (backwards_fd) {
            finite_difference_orders = {4};
        }
        else {
            finite_difference_orders = {2, 4};
        }

        for (size_t finite_difference_order : finite_difference_orders) {

            /* In this example we want to examine the convergence behavior under the assumption of exponential stay time 
    distributions. In this case, we can compare the solution of the IDE simulation with a corresponding ODE solution. */
            std::string save_dir;
            if (backwards_fd) {
                save_dir = fmt::format("../../simulation_results/"
                                       "detailed_init_exponential_late_t0/",
                                       ode_exponent, finite_difference_order);
            }
            else {
                save_dir = fmt::format("../../simulation_results/"
                                       "detailed_init_exponential_dt_ode=1e-{:.0f}_finite_diff={}_central_fd/",
                                       ode_exponent, finite_difference_order);
            }
            // Make folder if not existent yet.
            boost::filesystem::path dir(save_dir);
            boost::filesystem::create_directories(dir);

            auto result_ode = simulate_ode(ode_exponent, save_dir).value();

            // Do IDE simulations.
            std::vector<ScalarType> ide_exponents = {0, 1, 2, 3, 4};
            std::vector<size_t> gregory_orders    = {1, 2, 3};

            // std::vector<ScalarType> ide_exponents = {1};
            // std::vector<size_t> gregory_orders    = {1};

            for (size_t gregory_order : gregory_orders) {
                std::cout << std::endl;
                std::cout << "Gregory order: " << gregory_order << std::endl;
                mio::IOResult<void> result_ide = simulate_ide(ide_exponents, gregory_order, finite_difference_order,
                                                              save_dir, result_ode, backwards_fd);
            }
        }
    }
}
