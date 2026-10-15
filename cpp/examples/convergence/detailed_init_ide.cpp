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
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
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

/**
 * @brief Read a TimeSeries from a CSV file as written by TimeSeries::export_csv, i.e. with a header line and one row per
 * time point of the form "time,value_1,...,value_n".
 *
 * @param[in] filename Path to the CSV file.
 * @param[in] num_elements Number of values per time point (without the time column).
 * @return The TimeSeries containing the values of the CSV file.
 */
mio::IOResult<mio::TimeSeries<ScalarType>> read_timeseries_from_csv(const std::string& filename, size_t num_elements)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        return mio::failure(mio::StatusCode::FileNotFound, "Could not open file " + filename);
    }

    mio::TimeSeries<ScalarType> timeseries(num_elements);
    mio::TimeSeries<ScalarType>::Vector values(num_elements);

    std::string line;
    // Skip header.
    std::getline(file, line);
    while (std::getline(file, line)) {
        if (line.empty()) {
            continue;
        }
        std::stringstream line_stream(line);
        std::string entry;

        std::getline(line_stream, entry, ',');
        ScalarType time = std::stod(entry);
        for (size_t i = 0; i < num_elements; i++) {
            if (!std::getline(line_stream, entry, ',')) {
                return mio::failure(mio::StatusCode::InvalidFileFormat,
                                    "Too few columns in line \"" + line + "\" of file " + filename);
            }
            values[i] = std::stod(entry);
        }
        timeseries.add_time_point(time, values);
    }

    return mio::success(timeseries);
}

/**
 * @brief Get the index of the time point in timeseries that corresponds to time. The time points of timeseries are
 * assumed to be equidistant.
 *
 * Returns a failure if timeseries does not contain a time point that matches time (up to half of the step size of
 * timeseries), e.g. because time is outside of the time range of timeseries or the step size of timeseries is too
 * coarse.
 */
mio::IOResult<Eigen::Index> get_index_of_time(const mio::TimeSeries<ScalarType>& timeseries, ScalarType time)
{
    ScalarType t0 = timeseries.get_time(0);
    ScalarType dt = timeseries.get_time(1) - t0;

    Eigen::Index index = (Eigen::Index)std::round((time - t0) / dt);
    if (index < 0 || index >= timeseries.get_num_time_points() ||
        std::fabs(timeseries.get_time(index) - time) > dt / 2.) {
        return mio::failure(mio::StatusCode::OutOfRange,
                            fmt::format("Time point {} is not contained in the given initial values (t = {} ... {} "
                                        "with dt = {}).",
                                        time, t0, timeseries.get_last_time(), dt));
    }
    return mio::success(index);
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

    auto results = {compressed_compartments, compressed_flows};
    return mio::success(results);
}

mio::IOResult<void> simulate_ide(std::vector<ScalarType> ide_exponents, size_t gregory_order,
                                 size_t finite_difference_order, ScalarType t_init_window, ScalarType t0_ide,
                                 ScalarType tmax, ScalarType TimeInfected, ScalarType cont_freq,
                                 std::string save_dir = "", bool kahan = true,
                                 mio::TimeSeries<ScalarType> compartments_groundtruth =
                                     mio::TimeSeries<ScalarType>((size_t)mio::isir::InfectionState::Count),
                                 bool more_precise_s_deriv = false, bool forward_fd = false)
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

            // Initialize time points before t0_ide based on groundtruth. The values are taken from the groundtruth at
            // the time points t_init, t_init + dt_ide, ..., t0_ide, i.e. with the step size of the IDE model. This
            // requires that dt_ide is a multiple of the step size of the groundtruth.
            Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));

            ScalarType t_init = t0_ide - t_init_window;
            ScalarType t0_ode = compartments_groundtruth.get_time(0);
            cutoff_window     = t_init - t0_ode;

            ScalarType start_time = more_precise_s_deriv ? t0_ode : t_init;

            // Number of IDE time steps from start_time to t0_ide.
            size_t num_steps = (size_t)std::round((t0_ide - start_time) * div_dt_ide);
            for (size_t step = 0; step <= num_steps; step++) {
                // Compute time directly from step to avoid accumulating rounding errors.
                ScalarType time = start_time + step * dt_ide;
                BOOST_OUTCOME_TRY(auto&& index, get_index_of_time(compartments_groundtruth, time));
                vec_init = compartments_groundtruth.get_value(index);
                init_populations.add_time_point(time, vec_init);
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

        sim.advance(tmax, kahan, more_precise_s_deriv, forward_fd, cutoff_window);

        // sim.advance_S_deriv_analytical(tmax);

        if (!save_dir.empty()) {
            // Save compartments.
            mio::TimeSeries<ScalarType> compartments = sim.get_result();
            mio::TimeSeries<ScalarType> flows        = sim.get_flows();

            auto result = compartments.export_csv(fmt::format("{}/ide_result.csv", save_dir), {}, ',', 12);

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

    // CSV file with the values that are used to initialize the IDE model. The IDE model is initialized with the values
    // at the time points t_init, t_init + dt_ide, ..., t0_ide. Therefore, the CSV file has to cover this time range and
    // dt_ide has to be a multiple of the step size used in the CSV file.
    std::string init_values_file =
        "./simulation_results/2026-10-12/"
        "Init_values_dtode=1e-6_t0ode=0_timeinf=2_contfreq=0.73_kahan=false_buffer=false_forwardfd=true/"
        "detailed_init_exponential_t0ide=48_tmax=50_finite_diff=4/tinit=48/ide_result.csv";

    std::vector<ScalarType> time_infected_values = {2.};

    ScalarType t0_ode                    = 0.;
    ScalarType t0_ide                    = 50.;
    std::vector<ScalarType> init_windows = {2.};
    std::vector<ScalarType> tmax_values  = {t0_ide + 100.};

    bool kahan                = false;
    bool more_precise_s_deriv = false;
    bool forward_fd           = true;

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
                fmt::format("./simulation_results/2026-10-13/"
                            "init_ide_dtode=1e-{}_t0ode={}_timeinf={}_contfreq={}_kahan={}_buffer={}_forwardfd={}/"
                            "detailed_init_exponential_t0ide={}_tmax={}_finite_diff={}/",
                            ode_exponent, t0_ode, time_infected, cont_freq, kahan, more_precise_s_deriv, forward_fd,
                            t0_ide, tmax, finite_difference_order);

            // Make folder if not existent yet.
            std::filesystem::path dir(save_dir);
            std::filesystem::create_directories(dir);

            // ScalarType saving_exponent = *std::max_element(ide_exponents.begin(), ide_exponents.end());
            ScalarType saving_exponent = 3.;
            auto result_ode =
                simulate_ode(ode_exponent, t0_ode, tmax, time_infected, cont_freq, save_dir, saving_exponent).value();

            auto compartments_init_result =
                read_timeseries_from_csv(init_values_file, (size_t)mio::isir::InfectionState::Count);
            if (!compartments_init_result) {
                std::cout << compartments_init_result.error().formatted_message() << std::endl;
                return -1;
            }
            auto compartments_init = compartments_init_result.value();

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
                    mio::IOResult<void> result_ide = simulate_ide(
                        ide_exponents, gregory_order, finite_difference_order, init_window, t0_ide, tmax, time_infected,
                        cont_freq, save_dir_ide, kahan, compartments_init, more_precise_s_deriv, forward_fd);
                    if (!result_ide) {
                        std::cout << result_ide.error().formatted_message() << std::endl;
                        return -1;
                    }
                }
            }
        }
    }
}
