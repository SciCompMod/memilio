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
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "ode_sir/model.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/utils/time_series.h"
#include "memilio/io/result_io.h"
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <vector>

using namespace mio;
namespace params
{
size_t num_agegroups = 1;

ScalarType TransmissionProbabilityOnContact = 0.5;
ScalarType RiskOfInfectionFromSymptomatic   = 1.;
ScalarType Seasonality                      = 0.;

ScalarType cont_freq    = 1.;
ScalarType damping      = 0.1;
ScalarType scaling_time = 25.;

ScalarType S0               = 999000.;
ScalarType I0               = 1000.;
ScalarType R0               = 0.;
ScalarType total_population = S0 + I0 + R0;
} // namespace params

mio::UncertainContactMatrix<ScalarType> scale_contact_matrix(ScalarType damping, ScalarType scaling_time)
{
    using namespace params;

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup<ScalarType>(1, 1);
    if (damping <= 1.) {
        // Perform simulation with a decrease in contacts.
        contact_matrix[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixXd::Constant(1, 1, cont_freq));
        contact_matrix[0].add_damping(0., mio::SimulationTime<ScalarType>(scaling_time));
        contact_matrix[0].add_damping(damping, mio::SimulationTime<ScalarType>(scaling_time + 0.00001));
    }
    // else {
    //     // Perform simulation with an increase in contacts.
    //     contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, damping * cont_freq));
    //     contact_matrix[0].add_damping(1 - 1. / damping, mio::SimulationTime(-100.));
    //     contact_matrix[0].add_damping(1 - 1. / damping, mio::SimulationTime(scaling_time));
    //     contact_matrix[0].add_damping(0., mio::SimulationTime(scaling_time + 0.1));
    // }

    return mio::UncertainContactMatrix(contact_matrix);
}

mio::IOResult<mio::TimeSeries<ScalarType>> simulate_ode(ScalarType ode_exponent, ScalarType t0_ode, ScalarType tmax,
                                                        int TimeInfected, std::string save_dir = "")
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

    // mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    // contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));

    mio::UncertainContactMatrix<ScalarType> contact_matrix         = scale_contact_matrix(damping, scaling_time);
    model.parameters.get<mio::osir::ContactPatterns<ScalarType>>() = mio::UncertainContactMatrix(contact_matrix);

    model.check_constraints();

    std::unique_ptr<mio::OdeIntegratorCore<ScalarType>> integrator =
        std::make_unique<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_fehlberg78>>();

    auto sir = simulate<ScalarType, mio::osir::Model<ScalarType>>(t0_ode, tmax, dt_ode, model, std::move(integrator));

    if (!save_dir.empty()) {
        // Save compartments.
        mio::TimeSeries<ScalarType> compartments = sir;
        auto result                              = compartments.export_csv("ode_result.csv");
        auto save_result_status_ode =
            mio::save_result({compartments}, {0}, num_agegroups, save_dir + "result_ode" + ".h5");

        if (!save_result_status_ode) {
            return mio::failure(mio::StatusCode::InvalidValue,
                                "Error occured while saving the ODE simulation results.");
        }
    }
    return mio::success(sir);
}

mio::IOResult<void> simulate_ide(ScalarType ide_exponent, size_t gregory_order, size_t finite_difference_order,
                                 ScalarType t0_ode, ScalarType t0_ide, ScalarType tmax, ScalarType TimeInfected,
                                 std::string save_dir = "",
                                 mio::TimeSeries<ScalarType> result_groundtruth =
                                     mio::TimeSeries<ScalarType>((size_t)mio::isir::InfectionState::Count))
{
    using namespace params;
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType dt_ide = pow(10, -ide_exponent);
    std::cout << "Simulation with dt=" << dt_ide << std::endl;

    mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);

    if (result_groundtruth.get_num_time_points() == 0) {
        std::cout << "No groundtruth was given.\n";
    }
    else {
        std::cout << "Initializing with given groundtruth.\n";

        // Initialize time points before t0_ide based on groundtruth.
        ScalarType dt_groundtruth       = result_groundtruth.get_time(1) - result_groundtruth.get_time(0);
        size_t groundtruth_index_factor = size_t(dt_ide / dt_groundtruth);

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

    mio::StateAgeFunctionWrapper dist(exp);
    std::vector<mio::StateAgeFunctionWrapper<ScalarType>> vec_dist((size_t)mio::isir::InfectionTransition::Count, dist);
    model.parameters.get<mio::isir::TransitionDistributions>() = vec_dist;

    mio::ConstantFunction transmissiononcontact(TransmissionProbabilityOnContact);
    mio::StateAgeFunctionWrapper transmissiononcontact_wrapper(transmissiononcontact);
    model.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = transmissiononcontact_wrapper;

    mio::ConstantFunction riskofinfection(RiskOfInfectionFromSymptomatic);
    mio::StateAgeFunctionWrapper riskofinfection_wrapper(riskofinfection);
    model.parameters.get<mio::isir::RiskOfInfectionFromSymptomatic>() = riskofinfection_wrapper;

    // mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    // contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    mio::UncertainContactMatrix<ScalarType> contact_matrix = scale_contact_matrix(damping, scaling_time);
    model.parameters.get<mio::isir::ContactPatterns>()     = mio::UncertainContactMatrix(contact_matrix);

    std::cout << "support max: " << model.compute_calctime(dt_ide, 1e-8) << std::endl;

    // Carry out simulation.
    mio::isir::SimulationMessinaExtendedDetailedInit sim(model, dt_ide);
    sim.advance(tmax);

    if (!save_dir.empty()) {
        // Save compartments.
        mio::TimeSeries<ScalarType> compartments = sim.get_result();
        auto save_result_status_ide =
            mio::save_result({compartments}, {0}, num_agegroups, save_dir + "result_ide" + ".h5");

        if (!save_result_status_ide) {
            return mio::failure(mio::StatusCode::InvalidValue,
                                "Error occured while saving the IDE simulation results.");
        }
    }

    return mio::success();
}

int main()
{
    /* In this example we want to examine the convergence behavior under the assumption of exponential stay time 
    distributions. In this case, we can compare the solution of the IDE simulation with a corresponding ODE solution. */

    using namespace params;

    // Compute groundtruth with ODE model.
    ScalarType ode_exponent = 4;

    ScalarType time_infected = 4.;

    ScalarType t0_ode = 0.;
    ScalarType t0_ide = 0.;
    ScalarType tmax   = 30.;

    size_t finite_difference_order = 4;

    ScalarType ide_exponent = 1.;
    size_t gregory_order    = 3;

    std::string save_dir = fmt::format("../../simulation_results/2026-02-10/test/");

    // Make folder if not existent yet.
    boost::filesystem::path dir(save_dir);
    boost::filesystem::create_directories(dir);

    auto result_ode = simulate_ode(ode_exponent, t0_ode, tmax, time_infected, save_dir).value();

    // Do IDE simulations.
    std::cout << std::endl;
    std::cout << "Gregory order: " << gregory_order << std::endl;
    mio::IOResult<void> result_ide = simulate_ide(ide_exponent, gregory_order, finite_difference_order, t0_ode, t0_ide,
                                                  tmax, time_infected, save_dir, result_ode);
}