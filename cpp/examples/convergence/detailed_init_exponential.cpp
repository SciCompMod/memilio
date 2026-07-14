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
#include "memilio/compartments/simulation.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "ode_secir/model.h"
#include "ode_sir/model.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/utils/time_series.h"
#include "memilio/io/result_io.h"
#include <Eigen/src/Core/util/Meta.h>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <vector>

using namespace mio;
namespace params
{
size_t num_agegroups = 1;

ScalarType TransmissionProbabilityOnContact = 0.8;
ScalarType RiskOfInfectionFromSymptomatic   = 1.;
ScalarType Seasonality                      = 0.;

ScalarType cont_freq = 0.7;

ScalarType totalpop_reduction = 1.;
ScalarType total_population   = 1e7 / totalpop_reduction;
ScalarType I0                 = 1000. / totalpop_reduction;
ScalarType R0                 = 0. / totalpop_reduction;
ScalarType S0                 = total_population - I0 - R0;

} // namespace params

mio::TimeSeries<ScalarType> compress_timeseries(const mio::TimeSeries<ScalarType>& simulation_result,
                                                ScalarType saving_dt_exponent)
{
    mio::TimeSeries<ScalarType> removed(simulation_result.get_num_elements());
    ScalarType dt_original = simulation_result.get_time(1) - simulation_result.get_time(0);
    ScalarType time        = simulation_result.get_time(0); // =0
    ScalarType saving_dt   = std::pow(10, -saving_dt_exponent);
    for (int i = 0; i < simulation_result.get_num_time_points(); i++) {
        if (std::fabs(simulation_result.get_time(i) - time) < dt_original / 2.) {
            removed.add_time_point(simulation_result.get_time(i), simulation_result[i]);
            // std::cout << time << std::endl;
            time += saving_dt;
        }
    }
    return removed;
}

mio::IOResult<std::vector<mio::TimeSeries<ScalarType>>> simulate_ode(ScalarType ode_exponent, ScalarType t0_ode,
                                                                     ScalarType tmax, int TimeInfected,
                                                                     std::string save_dir       = "",
                                                                     ScalarType saving_exponent = 0.)
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

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup<ScalarType>(1, 1);
    contact_matrix[0]                      = mio::ContactMatrix<ScalarType>(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    // mio::UncertainContactMatrix<ScalarType> contact_matrix         = scale_contact_matrix(scaling_factor_contacts);
    model.parameters.get<mio::osir::ContactPatterns<ScalarType>>() = mio::UncertainContactMatrix(contact_matrix);

    model.check_constraints();

    std::unique_ptr<mio::OdeIntegratorCore<ScalarType>> integrator =
        std::make_unique<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_fehlberg78>>();

    auto sim = mio::FlowSimulation<ScalarType, mio::osir::Model<ScalarType>>(model, t0_ode, dt_ode);
    sim.set_integrator_core(std::move(integrator));
    sim.set_last_step_tolerance(1e-8);

    sim.advance(tmax);
    auto compartments = sim.get_result();
    auto flows        = sim.get_flows();

    std::cout << "Num tps ODE: " << compartments.get_num_time_points() << std::endl;

    mio::TimeSeries<ScalarType> compressed_compartments =
        mio::TimeSeries<ScalarType>((Eigen::Index)mio::isir::InfectionState::Susceptible);

    if (!save_dir.empty()) {
        // Save compartments.
        compressed_compartments = compress_timeseries(compartments, saving_exponent);

        std::cout << "Num tps ODE compressed: " << compressed_compartments.get_num_time_points() << std::endl;

        auto result = compressed_compartments.export_csv(fmt::format("{}/ode_result_compressed.csv", save_dir));
        auto save_result_status_ode =
            mio::save_result({compressed_compartments}, {0}, num_agegroups,
                             save_dir + "result_ode_dt=1e-" + fmt::format("{:.0f}", ode_exponent) + "_savedt=1e-" +
                                 fmt::format("{:.0f}", saving_exponent) + ".h5");

        if (!save_result_status_ode) {
            return mio::failure(mio::StatusCode::InvalidValue,
                                "Error occured while saving the ODE simulation results.");
        }
    }

    auto results = {compressed_compartments, flows};
    return mio::success(results);
}

mio::IOResult<void> simulate_ide(std::vector<ScalarType> ide_exponents, ScalarType ode_exponent, size_t gregory_order,
                                 size_t finite_difference_order, ScalarType t_init_window, ScalarType t0_ide,
                                 ScalarType tmax, ScalarType TimeInfected, std::string save_dir = "",
                                 mio::TimeSeries<ScalarType> compartments_groundtruth =
                                     mio::TimeSeries<ScalarType>((size_t)mio::isir::InfectionState::Count),
                                 mio::TimeSeries<ScalarType> flows_groundtruth =
                                     mio::TimeSeries<ScalarType>((size_t)mio::isir::InfectionTransition::Count))
{
    using namespace params;
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    for (ScalarType ide_exponent : ide_exponents) {

        ScalarType dt_ide     = pow(10, -ide_exponent);
        ScalarType div_dt_ide = pow(10, ide_exponent);
        std::cout << "Simulation with dt=" << dt_ide << std::endl;

        mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);
        mio::TimeSeries<ScalarType> init_flows_ts((size_t)mio::isir::InfectionTransition::Count);

        if (compartments_groundtruth.get_num_time_points() == 0) {
            std::cout << "No groundtruth was given.\n";
        }
        else {
            std::cout << "Initializing with given groundtruth for compartments.\n";

            // Initialize time points before t0_ide based on groundtruth.
            ScalarType div_dt_groundtruth = std::pow(10, ode_exponent);
            // Compute scaling of ode_exponent/ide_exponent or dt_ide/dt_ode.
            ScalarType groundtruth_index_factor = std::pow(10, ode_exponent - ide_exponent);
            std::cout << "groundtruth_index_factor: " << groundtruth_index_factor << std::endl;

            Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));
            Vec vec_init_flows(Vec::Constant((size_t)mio::isir::InfectionTransition::Count, 0.));

            std::vector<size_t> compartments = {(size_t)mio::isir::InfectionState::Susceptible,
                                                (size_t)mio::isir::InfectionState::Infected,
                                                (size_t)mio::isir::InfectionState::Recovered};

            std::vector<size_t> flows = {(size_t)mio::isir::InfectionTransition::SusceptibleToInfected,
                                         (size_t)mio::isir::InfectionTransition::InfectedToRecovered};

            // ScalarType t0_ode = compartments_groundtruth.get_time(0);
            ScalarType t_init = t0_ide - t_init_window;
            int start_index   = std::round(t_init * div_dt_groundtruth);

            // Add values to init_populations.
            for (size_t compartment : compartments) {
                vec_init[compartment] = compartments_groundtruth.get_value(start_index)[compartment];
            }
            std::cout << "t init index: " << start_index << std::endl;

            init_populations.add_time_point(t_init, vec_init);

            std::cout << "init_pop t_init: " << init_populations.get_last_time() << std::endl;

            std::cout << "diff S at t_init: "
                      << compartments_groundtruth.get_value(
                             std::round(t_init * div_dt_groundtruth))[(size_t)mio::isir::InfectionState::Susceptible] -
                             init_populations.get_last_value()[(size_t)mio::isir::InfectionState::Susceptible]
                      << std::endl;

            while (init_populations.get_last_time() < t0_ide - 1e-10) {
                start_index += groundtruth_index_factor;
                for (size_t compartment : compartments) {
                    vec_init[compartment] = compartments_groundtruth.get_value(std::round(start_index))[compartment];
                }

                init_populations.add_time_point(init_populations.get_last_time() + dt_ide, vec_init);
            }

            if (flows_groundtruth.get_num_time_points() > 0) {
                std::cout << "Initializing with given groundtruth for flows.\n";
                // Add values to init_flows_ts.
                vec_init_flows[(size_t)mio::isir::InfectionTransition::SusceptibleToInfected] =
                    compartments_groundtruth.get_value(0)[(size_t)mio::isir::InfectionState::Infected];
                vec_init_flows[(size_t)mio::isir::InfectionTransition::InfectedToRecovered] =
                    compartments_groundtruth.get_value(0)[(size_t)mio::isir::InfectionState::Recovered];

                init_flows_ts.add_time_point(t_init, vec_init_flows);

                while (init_flows_ts.get_last_time() < t_init - 1e-10) {
                    for (size_t flow : flows) {
                        vec_init_flows[flow] =
                            (flows_groundtruth.get_value(
                                 int(std::round(t_init * div_dt_groundtruth) +
                                     init_flows_ts.get_num_time_points() * groundtruth_index_factor))[flow] -
                             flows_groundtruth.get_value(
                                 int(std::round(t_init * div_dt_groundtruth) +
                                     (init_flows_ts.get_num_time_points() - 1) * groundtruth_index_factor))[flow]) /
                            dt_ide;
                    }
                    init_flows_ts.add_time_point(init_flows_ts.get_last_time() + dt_ide, vec_init_flows);
                }
            }
        }

        // Initialize model.
        mio::isir::ModelMessinaExtendedDetailedInit model(std::move(init_populations), total_population, gregory_order,
                                                          finite_difference_order, std::move(init_flows_ts));

        mio::ExponentialSurvivalFunction exp(1. / TimeInfected);

        mio::StateAgeFunctionWrapper dist(exp);
        std::vector<mio::StateAgeFunctionWrapper<ScalarType>> vec_dist((size_t)mio::isir::InfectionTransition::Count,
                                                                       dist);
        model.parameters.get<mio::isir::TransitionDistributions>() = vec_dist;

        mio::ConstantFunction transmissiononcontact(TransmissionProbabilityOnContact);
        mio::StateAgeFunctionWrapper transmissiononcontact_wrapper(transmissiononcontact);
        model.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = transmissiononcontact_wrapper;

        mio::ConstantFunction riskofinfection(RiskOfInfectionFromSymptomatic);
        mio::StateAgeFunctionWrapper riskofinfection_wrapper(riskofinfection);
        model.parameters.get<mio::isir::RiskOfInfectionFromSymptomatic>() = riskofinfection_wrapper;

        mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup<ScalarType>(1, 1);
        contact_matrix[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixXd::Constant(1, 1, cont_freq));
        // mio::UncertainContactMatrix<ScalarType> contact_matrix = scale_contact_matrix(scaling_factor_contacts);
        model.parameters.get<mio::isir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

        std::cout << "support max: " << model.compute_calctime(dt_ide, 1e-7) << std::endl;

        // Carry out simulation.
        mio::isir::SimulationMessinaExtendedDetailedInit sim(model, dt_ide, div_dt_ide);
        // size_t fd_order_contacts = 1;
        sim.advance(tmax);
        // sim.advance_S_deriv_fixedpoint(tmax);

        if (!save_dir.empty()) {
            // Save compartments.
            mio::TimeSeries<ScalarType> compartments = sim.get_result();
            mio::TimeSeries<ScalarType> flows        = sim.get_flows();

            auto result = compartments.export_csv(fmt::format("{}/ide_result.csv", save_dir));

            auto save_result_status_ide =
                mio::save_result({compartments}, {0}, num_agegroups,
                                 save_dir + "result_ide_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
                                     "_gregoryorder=" + fmt::format("{}", gregory_order) + ".h5");
            auto save_result_status_ide_flows =
                mio::save_result({flows}, {0}, num_agegroups,
                                 save_dir + "result_ide_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
                                     "_gregoryorder=" + fmt::format("{}", gregory_order) + "_flows.h5");

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

    // Compute groundtruth with ODE model.
    ScalarType ode_exponent = 7.;

    std::vector<ScalarType> time_infected_values = {2};
    // Support max with tol = 1e-8:
    // T_I=1: support max = 18.43
    // T_I=2: support max = 36.85
    // T_I=3: support max = 55.27
    // T_I=4: support max = 73.69

    ScalarType t0_ode                      = 0.;
    ScalarType t0_ide                      = 40.;
    std::vector<ScalarType> t_init_windows = {30.};
    std::vector<ScalarType> tmax_values    = {t0_ide + 5.};

    std::vector<size_t> finite_difference_orders = {4};

    std::vector<ScalarType> ide_exponents = {0, 1, 2, 3};
    std::vector<size_t> gregory_orders    = {1, 2, 3};

    std::vector<std::vector<ScalarType>> timeinf_tmax_values;

    for (ScalarType time_infected : time_infected_values) {
        for (ScalarType tmax : tmax_values) {
            std::vector<ScalarType> value_tuple = {time_infected, tmax};
            timeinf_tmax_values.push_back(value_tuple);
        }
    }

    for (std::vector<ScalarType> value_tuple : timeinf_tmax_values) {

        ScalarType time_infected = value_tuple[0];
        ScalarType tmax          = value_tuple[1];

        for (size_t finite_difference_order : finite_difference_orders) {
            std::cout << "FD order: " << finite_difference_order << std::endl;

            std::string save_dir =
                fmt::format("../../simulation_results/2026-07-13/more_kahan_dtode=1e-{}_t0ode={}_timeinf={}/"
                            "detailed_init_exponential_t0ide={}_tmax={}_finite_diff={}/",
                            ode_exponent, t0_ode, time_infected, t0_ide, tmax, finite_difference_order);

            // Make folder if not existent yet.
            std::filesystem::path dir(save_dir);
            std::filesystem::create_directories(dir);

            ScalarType saving_exponent = *std::max_element(ide_exponents.begin(), ide_exponents.end());
            // ScalarType saving_exponent = ode_exponent;
            auto result_ode =
                simulate_ode(ode_exponent, t0_ode, tmax, time_infected, save_dir, saving_exponent).value();

            auto compartments_ode = result_ode[0];
            // auto flows_ode        = result_ode[1];

            for (ScalarType t_init_window : t_init_windows) {
                ScalarType t_init = t0_ide - t_init_window;
                std::cout << "t_init = " << t_init << std::endl;

                std::string save_dir_ide = fmt::format("{}/tinit={}/", save_dir, t_init);
                // Make folder if not existent yet.
                std::filesystem::path dir_ide(save_dir_ide);
                std::filesystem::create_directories(dir_ide);

                // Do IDE simulations.
                for (size_t gregory_order : gregory_orders) {
                    std::cout << std::endl;
                    std::cout << "Gregory order: " << gregory_order << std::endl;
                    mio::IOResult<void> result_ide =
                        simulate_ide(ide_exponents, saving_exponent, gregory_order, finite_difference_order,
                                     t_init_window, t0_ide, tmax, time_infected, save_dir_ide, compartments_ode);
                }
            }
        }
    }
}
