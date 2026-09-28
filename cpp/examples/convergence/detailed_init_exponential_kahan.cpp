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

ScalarType kahan_sum(const std::vector<ScalarType>& v)
{
    long double sum = 0.0L;
    long double c   = 0.0L;
    for (auto x : v) {
        long double y = (long double)x - c;
        long double t = sum + y;
        c             = (t - sum) - y;
        sum           = t;
    }
    return static_cast<ScalarType>(sum);
}

// First order backward finite difference approximation of S'(t_index) on the (raw, dt_ode-spaced) groundtruth
// Susceptible values, mirroring ModelMessinaExtendedDetailedInit::compute_S_deriv's first order stencil. Same
// cancellation concern as the fourth order stencil below, so the two terms are also Kahan-summed.
ScalarType s_deriv_fd1_kahan(const mio::TimeSeries<ScalarType>& groundtruth, int index, int stride, ScalarType div_dt)
{
    std::vector<ScalarType> terms = {
        groundtruth.get_value(index)[(size_t)mio::isir::InfectionState::Susceptible],
        -groundtruth.get_value(index - stride)[(size_t)mio::isir::InfectionState::Susceptible]};
    return -kahan_sum(terms) * div_dt;
}

// Fourth order backward finite difference approximation of S'(t_index) on the (raw, dt_ode-spaced) groundtruth
// Susceptible values, mirroring ModelMessinaExtendedDetailedInit::compute_S_deriv's fourth order stencil.
// The stencil differences values of order total_population (~1e7) over a very small dt_ode, so the terms nearly
// cancel; summing them in double precision loses most significant digits before the division by div_dt/12
// amplifies whatever rounding error remains. Accumulating with Kahan (compensated) summation in long double,
// as mio::isir::kahan_sum does for the analogous cancellation-prone sums in model.cpp, keeps that rounding error
// far below the O(dt_ode^4) truncation error the stencil is supposed to have.
ScalarType s_deriv_fd4_kahan(const mio::TimeSeries<ScalarType>& groundtruth, int index, int stride, ScalarType div_dt)
{
    std::vector<ScalarType> terms = {
        25 * groundtruth.get_value(index)[(size_t)mio::isir::InfectionState::Susceptible],
        -48 * groundtruth.get_value(index - stride)[(size_t)mio::isir::InfectionState::Susceptible],
        36 * groundtruth.get_value(index - 2 * stride)[(size_t)mio::isir::InfectionState::Susceptible],
        -16 * groundtruth.get_value(index - 3 * stride)[(size_t)mio::isir::InfectionState::Susceptible],
        3 * groundtruth.get_value(index - 4 * stride)[(size_t)mio::isir::InfectionState::Susceptible]};
    return -kahan_sum(terms) * (div_dt / 12.);
}

// Dispatches to the backward finite difference stencil of the requested order (1 or 4) for S'(t_index). The stencil
// points are spaced `stride` groundtruth indices apart, i.e. the effective step size is stride * dt_ode = 1/div_dt.
ScalarType s_deriv_fd_kahan(const mio::TimeSeries<ScalarType>& groundtruth, int index, int stride, ScalarType div_dt,
                            size_t order)
{
    switch (order) {
    case 1:
        return s_deriv_fd1_kahan(groundtruth, index, stride, div_dt);
    case 4:
        return s_deriv_fd4_kahan(groundtruth, index, stride, div_dt);
    default:
        throw std::invalid_argument("s_deriv_fd_kahan: unsupported finite difference order (must be 1 or 4).");
    }
}

mio::IOResult<std::vector<mio::TimeSeries<ScalarType>>> simulate_ode(ScalarType ode_exponent, ScalarType t0_ode,
                                                                     ScalarType tmax, ScalarType TimeInfected,
                                                                     ScalarType cont_freq, std::string save_dir = "",
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

    mio::TimeSeries<ScalarType> compressed_compartments = compress_timeseries(compartments, saving_exponent);
    mio::TimeSeries<ScalarType> compressed_flows        = compress_timeseries(flows, saving_exponent);

    std::cout << "Num tps ODE compressed: " << compressed_compartments.get_num_time_points() << std::endl;

    if (!save_dir.empty()) {
        // Save compartments.
        // auto result = compressed_compartments.export_csv(fmt::format("{}/ode_result_compressed.csv", save_dir));
        auto save_result_status_ode =
            mio::save_result({compressed_compartments}, {0}, num_agegroups,
                             save_dir + "result_ode_dt=1e-" + fmt::format("{:.0f}", ode_exponent) + "_savedt=1e-" +
                                 fmt::format("{:.0f}", saving_exponent) + ".h5");

        auto save_result_status_ode_flows =
            mio::save_result({compressed_flows}, {0}, num_agegroups,
                             save_dir + "result_ode_dt=1e-" + fmt::format("{:.0f}", ode_exponent) + "_savedt=1e-" +
                                 fmt::format("{:.0f}", saving_exponent) + "_flows.h5");

        if (!save_result_status_ode) {
            return mio::failure(mio::StatusCode::InvalidValue,
                                "Error occured while saving the ODE simulation results.");
        }
    }

    auto results = {compartments, flows};
    return mio::success(results);
}

mio::IOResult<void> simulate_ide(std::vector<ScalarType> ide_exponents, ScalarType ode_exponent, size_t gregory_order,
                                 size_t finite_difference_order, ScalarType t_init_window, ScalarType t0_ide,
                                 ScalarType tmax, ScalarType TimeInfected, ScalarType cont_freq,
                                 std::string save_dir = "", bool kahan = true,
                                 mio::TimeSeries<ScalarType> compartments_groundtruth =
                                     mio::TimeSeries<ScalarType>((size_t)mio::isir::InfectionState::Count),
                                 bool more_precise_s_deriv = false, bool forward_fd = false,
                                 mio::TimeSeries<ScalarType> flows_groundtruth =
                                     mio::TimeSeries<ScalarType>((size_t)mio::isir::InfectionState::Count),
                                 size_t s_deriv_fd_order = 0, ScalarType flow_init_exponent = -1.)
{
    using namespace params;
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    for (ScalarType ide_exponent : ide_exponents) {

        ScalarType dt_ide     = pow(10, -ide_exponent);
        ScalarType div_dt_ide = pow(10, ide_exponent);
        std::cout << "Simulation with dt=" << dt_ide << std::endl;

        mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);
        mio::TimeSeries<ScalarType> init_flows_ts((size_t)mio::isir::InfectionTransition::Count);

        ScalarType cutoff_window = 0;

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

            ScalarType t_init = t0_ide - t_init_window;
            ScalarType t0_ode = compartments_groundtruth.get_time(0);
            cutoff_window     = t_init - t0_ode;

            int start_index       = 0;
            ScalarType start_time = 0;
            if (more_precise_s_deriv) {
                start_time = t0_ode;
            }
            else {
                start_time = t_init;
            }
            // compartments_groundtruth is indexed from its own start time t0_ode, not from absolute time 0, so the
            // index into it has to be computed relative to t0_ode.
            start_index = std::round((start_time - t0_ode) * div_dt_groundtruth);

            // Add values to init_populations.
            for (size_t compartment : compartments) {
                vec_init[compartment] = compartments_groundtruth.get_value(start_index)[compartment];
            }

            init_populations.add_time_point(start_time, vec_init);

            while (init_populations.get_last_time() < t0_ide - 1e-10) {
                std::round(start_index += groundtruth_index_factor);
                for (size_t compartment : compartments) {
                    vec_init[compartment] = compartments_groundtruth.get_value(std::round(start_index))[compartment];
                }

                init_populations.add_time_point(init_populations.get_last_time() + dt_ide, vec_init);
            }

            if (flows_groundtruth.get_num_time_points() > 0) {
                std::cout << "Initializing with given groundtruth for flows.\n";
                // Compute the flow rate at t_init as a backward difference of the cumulative groundtruth flows,
                // consistent with how every later point of init_flows_ts is derived below. (Using compartment
                // values here, as was done previously, mixes up populations and flow rates and leaves the point
                // at t_init inconsistent with the rest of the series.)
                // The backward differences use the step size dt_flow_init = 10^-flow_init_exponent (default: dt_ode),
                // which is independent of dt_ide. It has to be a multiple of dt_ode so that all stencil points
                // coincide with groundtruth time points. Note that very small steps amplify rounding errors.
                ScalarType used_flow_init_exponent = flow_init_exponent < 0. ? ode_exponent : flow_init_exponent;
                if (used_flow_init_exponent > ode_exponent) {
                    return mio::failure(mio::StatusCode::InvalidValue,
                                        "flow_init_exponent must not exceed ode_exponent, i.e. the flow "
                                        "initialization step must not be smaller than dt_ode.");
                }
                ScalarType div_dt_flow_init = std::pow(10, used_flow_init_exponent);
                int flow_init_stride        = int(std::round(std::pow(10, ode_exponent - used_flow_init_exponent)));
                std::cout << "Flow initialization with dt = 1e-" << used_flow_init_exponent
                          << " (stride in groundtruth indices: " << flow_init_stride << ")" << std::endl;
                int t_init_index = int(std::round(t_init * div_dt_groundtruth));
                for (size_t flow : flows) {
                    if (s_deriv_fd_order > 0 && flow == (size_t)mio::isir::InfectionTransition::SusceptibleToInfected) {
                        // Compute the S -> I flow as -S' via a backward finite difference scheme of the requested
                        // order, applied to the groundtruth Susceptible values with step dt_flow_init.
                        vec_init_flows[flow] = s_deriv_fd_kahan(compartments_groundtruth, t_init_index,
                                                                flow_init_stride, div_dt_flow_init, s_deriv_fd_order);
                    }
                    else {
                        vec_init_flows[flow] = (flows_groundtruth.get_value(t_init_index)[flow] -
                                                flows_groundtruth.get_value(t_init_index - flow_init_stride)[flow]) *
                                               div_dt_flow_init;
                    }
                }

                init_flows_ts.add_time_point(t_init, vec_init_flows);

                while (init_flows_ts.get_last_time() < t0_ide - 1e-10) {
                    int index = int(std::round(t_init * div_dt_groundtruth) +
                                    init_flows_ts.get_num_time_points() * groundtruth_index_factor);
                    for (size_t flow : flows) {
                        if (s_deriv_fd_order > 0 &&
                            flow == (size_t)mio::isir::InfectionTransition::SusceptibleToInfected) {
                            vec_init_flows[flow] = s_deriv_fd_kahan(compartments_groundtruth, index, flow_init_stride,
                                                                    div_dt_flow_init, s_deriv_fd_order);
                        }
                        else {
                            vec_init_flows[flow] = (flows_groundtruth.get_value(index)[flow] -
                                                    flows_groundtruth.get_value(index - flow_init_stride)[flow]) *
                                                   div_dt_flow_init;
                        }
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

        // std::cout << "support max: " << model.compute_calctime(dt_ide, 1e-7) << std::endl;

        // Carry out simulation.
        mio::isir::SimulationMessinaExtendedDetailedInit sim(model, dt_ide, div_dt_ide);
        // size_t fd_order_contacts = 1;

        sim.advance(tmax, kahan, more_precise_s_deriv, forward_fd, cutoff_window);

        // sim.advance_S_deriv_analytical(tmax);

        if (!save_dir.empty()) {
            // Save compartments.
            mio::TimeSeries<ScalarType> compartments = sim.get_result();
            mio::TimeSeries<ScalarType> flows        = sim.get_flows();

            // auto result = compartments.export_csv(fmt::format("{}/ide_result.csv", save_dir));

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
    ScalarType ode_exponent = 6.;

    std::vector<ScalarType> time_infected_values = {2.};

    ScalarType t0_ode                    = 0.;
    ScalarType t0_ide                    = 50.;
    std::vector<ScalarType> init_windows = {0., 10.};
    std::vector<ScalarType> tmax_values  = {t0_ide + 100.};

    bool kahan                = false;
    bool more_precise_s_deriv = false;
    bool forward_fd           = true;
    // Order of the backward finite difference scheme used to derive the initial S -> I flow from the groundtruth
    // Susceptible values; 0 disables this and falls back to differencing the groundtruth's cumulative flows.
    // size_t s_deriv_fd_order = 0;
    // Step size dt = 1e-flow_init_exponent used in the finite differences that initialize the flows. Must not exceed
    // ode_exponent (dt_ode is the finest possible choice); a negative value selects dt_ode.
    // ScalarType flow_init_exponent = 6.;

    std::vector<size_t> finite_difference_orders = {4};

    std::vector<ScalarType> ide_exponents = {3.};
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

        ScalarType cont_freq = 1.46 / time_infected;

        for (size_t finite_difference_order : finite_difference_orders) {
            std::cout << "FD order: " << finite_difference_order << std::endl;

            std::string save_dir =
                fmt::format("./simulation_results/2026-09-24/"
                            "forward_ta_dtode=1e-{}_t0ode={}_timeinf={}_contfreq={}_kahan={}_buffer={}_forwardfd={}/"
                            "detailed_init_exponential_t0ide={}_tmax={}_finite_diff={}/",
                            ode_exponent, t0_ode, time_infected, cont_freq, kahan, more_precise_s_deriv, forward_fd,
                            t0_ide, tmax, finite_difference_order);

            // Make folder if not existent yet.
            std::filesystem::path dir(save_dir);
            std::filesystem::create_directories(dir);

            // ScalarType saving_exponent = *std::max_element(ide_exponents.begin(), ide_exponents.end());
            ScalarType saving_exponent = 3.;
            // ScalarType saving_exponent = ode_exponent;
            auto result_ode =
                simulate_ode(ode_exponent, t0_ode, tmax, time_infected, cont_freq, save_dir, saving_exponent).value();

            auto compartments_ode = result_ode[0];
            auto flows_ode        = result_ode[1];

            for (ScalarType init_window : init_windows) {
                ScalarType t_init = t0_ide - init_window;
                std::cout << "t_init = " << t_init << std::endl;

                std::string save_dir_ide = fmt::format("{}/tinit={}/", save_dir, t_init);
                // Make folder if not existent yet.
                std::filesystem::path dir_ide(save_dir_ide);
                std::filesystem::create_directories(dir_ide);

                // Do IDE simulations.
                for (size_t gregory_order : gregory_orders) {
                    std::cout << std::endl;
                    std::cout << "Gregory order: " << gregory_order << std::endl;
                    // compartments_ode and flows_ode are raw (dt_ode-spaced) groundtruth series, so the
                    // groundtruth resolution passed to simulate_ide must be the true ode_exponent, not
                    // saving_exponent (which only applies to the compressed results written to disk).
                    mio::IOResult<void> result_ide =
                        simulate_ide(ide_exponents, ode_exponent, gregory_order, finite_difference_order, init_window,
                                     t0_ide, tmax, time_infected, cont_freq, save_dir_ide, kahan, compartments_ode,
                                     more_precise_s_deriv, forward_fd);
                    //flows_ode, s_deriv_fd_order, flow_init_exponent
                }
            }
        }
    }
}
