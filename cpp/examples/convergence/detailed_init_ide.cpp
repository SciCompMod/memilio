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
#include <algorithm>
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

/**
 * @brief Get the initial values at t0 from the given TimeSeries, e.g. from the compartments of the ODE solution.
 *
 * Using the ODE solution at t0 ensures that the IDE solution can be compared with the ODE solution at the same time
 * points also if the IDE simulation starts later than the ODE simulation.
 *
 * @param[in] values TimeSeries with equidistant time points that contains t0.
 * @param[in] t0 Start time of the IDE simulation.
 * @return TimeSeries with one time point t0 containing the initial values.
 */
mio::IOResult<mio::TimeSeries<ScalarType>> get_init_values_at_t0(const mio::TimeSeries<ScalarType>& values,
                                                                 ScalarType t0)
{
    BOOST_OUTCOME_TRY(auto&& index, get_index_of_time(values, t0));

    mio::TimeSeries<ScalarType> init_populations(values.get_num_elements());
    init_populations.add_time_point(t0, values.get_value(index));
    return mio::success(init_populations);
}

/**
 * @brief Get the initial values for the IDE model from the given TimeSeries.
 *
 * The values are taken from the given TimeSeries at the time points t_init, t_init + dt, ..., t0 with t_init = t0 -
 * t_init_window, i.e. with the step size of the IDE model. This requires that dt is a multiple of the step size of the
 * given TimeSeries. This can be used for both compartments and flows.
 *
 * @param[in] values TimeSeries with the values that are used for the initialization.
 * @param[in] t_init_window Length of the time window before t0 for which initial values are given.
 * @param[in] t0 Start time of the IDE simulation.
 * @param[in] dt Step size of the IDE simulation.
 * @param[in] div_dt Inverse of dt.
 * @return TimeSeries with the initial values for the IDE model.
 */
mio::IOResult<mio::TimeSeries<ScalarType>> get_init_values(const mio::TimeSeries<ScalarType>& values,
                                                           ScalarType t_init_window, ScalarType t0, ScalarType dt,
                                                           ScalarType div_dt)
{
    mio::TimeSeries<ScalarType> init_values(values.get_num_elements());

    if (values.get_num_time_points() == 0) {
        std::cout << "No initial values were given.\n";
        return mio::success(init_values);
    }

    ScalarType t_init = t0 - t_init_window;

    // Number of IDE time steps from t_init to t0.
    size_t num_steps = (size_t)std::round(t_init_window * div_dt);
    for (size_t step = 0; step <= num_steps; step++) {
        // Compute time directly from step to avoid accumulating rounding errors.
        ScalarType time = t_init + step * dt;
        BOOST_OUTCOME_TRY(auto&& index, get_index_of_time(values, time));
        init_values.add_time_point(time, values.get_value(index));
    }

    return mio::success(init_values);
}

/**
 * @brief Set up an IDE model with the given initial values and simulate it until tmax with step size dt.
 *
 * If init_flows is empty, the flows (i.e. -S') on the initial time window are computed by the simulation from
 * init_populations with finite differences. Otherwise, the given flows are used and have to be given at the same time
 * points as init_populations.
 *
 * @return Vector containing the compartments and the flows of the simulation result.
 */
std::vector<mio::TimeSeries<ScalarType>>
run_ide(mio::TimeSeries<ScalarType>&& init_populations, mio::TimeSeries<ScalarType>&& init_flows, ScalarType dt,
        ScalarType div_dt, size_t gregory_order, size_t finite_difference_order, ScalarType tmax,
        ScalarType TimeInfected, ScalarType cont_freq, bool kahan, bool forward_fd)
{
    using namespace params;

    // Initialize model.
    mio::isir::ModelMessinaExtendedDetailedInit model(std::move(init_populations), total_population, gregory_order,
                                                      finite_difference_order, std::move(init_flows));

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

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup<ScalarType>(1, 1);
    contact_matrix[0]                      = mio::ContactMatrix<ScalarType>(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    // mio::UncertainContactMatrix<ScalarType> contact_matrix = scale_contact_matrix(scaling_factor_contacts);
    model.parameters.get<mio::isir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    // std::cout << "support max: " << model.compute_calctime(dt, 1e-7) << std::endl;

    // Carry out simulation.
    mio::isir::SimulationMessinaExtendedDetailedInit sim(model, dt, div_dt);

    sim.advance(tmax, kahan, forward_fd);

    return {sim.get_result(), sim.get_flows()};
}

/**
 * @brief Simulate the IDE model from t0_ide until t_switch with the (smaller) step size dt_fine = 10^(-fine_exponent)
 * and save the resulting compartments and flows in CSV files.
 *
 * The simulation does not use an initial time window, i.e. it is only initialized with the values init_values_t0 at
 * t0_ide. The CSV files cover the time points t0_ide, ..., t_switch with step size dt_fine. The flows are given by
 * -S', where S' is computed with finite differences with step size dt_fine. The files can be read in with
 * read_timeseries_from_csv() to initialize a subsequent IDE simulation that starts at t_switch.
 *
 * @param[in] init_values_t0 TimeSeries with one time point t0_ide containing the initial values.
 * @param[in] csv_file Path to the CSV file for the compartments.
 * @param[in] csv_file_flows Path to the CSV file for the flows.
 */
mio::IOResult<void> simulate_ide_until_switch(ScalarType fine_exponent, ScalarType t_switch, size_t gregory_order,
                                              size_t finite_difference_order,
                                              const mio::TimeSeries<ScalarType>& init_values_t0,
                                              ScalarType TimeInfected, ScalarType cont_freq,
                                              const std::string& csv_file, const std::string& csv_file_flows,
                                              bool kahan, bool forward_fd)
{
    ScalarType dt_fine     = pow(10, -fine_exponent);
    ScalarType div_dt_fine = pow(10, fine_exponent);
    std::cout << "Simulation until t_switch=" << t_switch << " with dt_fine=" << dt_fine << std::endl;

    // Only use the values at t0_ide as initial values, i.e. no initial time window.
    std::vector<mio::TimeSeries<ScalarType>> result =
        run_ide(mio::TimeSeries<ScalarType>(init_values_t0),
                mio::TimeSeries<ScalarType>((size_t)mio::isir::InfectionTransition::Count), dt_fine, div_dt_fine,
                gregory_order, finite_difference_order, t_switch, TimeInfected, cont_freq, kahan, forward_fd);

    // Use a high precision so that the subsequent simulation is not affected by rounding of the saved values.
    BOOST_OUTCOME_TRY(result[0].export_csv(csv_file, {}, ',', 16));
    BOOST_OUTCOME_TRY(result[1].export_csv(csv_file_flows, {}, ',', 16));

    return mio::success();
}

/**
 * @brief Simulate the IDE model with the target step sizes dt = 10^(-ide_exponent) for all ide_exponents.
 *
 * For each target step size dt, the simulation is split into two parts. First, the IDE model is simulated on
 * [t0_ide, t_switch] with the smaller step size dt_fine = dt * 10^(-fine_exponent_offset), where
 * t_switch = t0_ide + num_steps_until_switch * dt and fine_exponent_offset is the entry of fine_exponent_offsets that
 * corresponds to ide_exponent. This simulation does not use an initial time window, i.e. it is only initialized with
 * the values of ode_compartments at t0_ide. The resulting compartments and flows are saved in CSV files. Afterwards,
 * the CSV files are read in and used to initialize the IDE model that is then simulated on [t_switch, tmax] with the
 * target step size dt. The initial time window of this simulation is the whole time interval [t0_ide, t_switch] of
 * the first simulation. If use_fine_flows, the flows (i.e. -S') computed with step
 * size dt_fine are used on the initial time window of the second simulation. Otherwise, they are computed from the
 * compartments with step size dt. Both parts use the same numerical scheme.
 * If num_steps_until_switch is zero, there is no simulation with step size dt_fine. Then, the IDE model is directly
 * simulated on [t0_ide, tmax] with step size dt and is only initialized with the values of ode_compartments at t0_ide.
 *
 * @param[in] ode_compartments Compartments of the ODE solution. The IDE model is initialized with its values at t0_ide
 *  so that the IDE solution can be compared with the ODE solution at the same time points.
 */
mio::IOResult<void> simulate_ide(std::vector<ScalarType> ide_exponents, std::vector<ScalarType> fine_exponent_offsets,
                                 size_t gregory_order, size_t finite_difference_order,
                                 const mio::TimeSeries<ScalarType>& ode_compartments, ScalarType t0_ide,
                                 ScalarType tmax, ScalarType TimeInfected, ScalarType cont_freq,
                                 size_t num_steps_until_switch, bool use_fine_flows, std::string save_dir = "",
                                 bool kahan = true, bool forward_fd = false)
{
    using namespace params;

    if (fine_exponent_offsets.size() != ide_exponents.size()) {
        return mio::failure(mio::StatusCode::InvalidValue,
                            "fine_exponent_offsets has to contain one entry for each entry of ide_exponents.");
    }

    BOOST_OUTCOME_TRY(auto&& init_values_t0, get_init_values_at_t0(ode_compartments, t0_ide));

    for (size_t i = 0; i < ide_exponents.size(); i++) {
        ScalarType ide_exponent         = ide_exponents[i];
        ScalarType fine_exponent_offset = fine_exponent_offsets[i];

        // dt has to be a multiple of dt_fine so that the values at the time points of the coarse grid are contained in
        // the result of the simulation with dt_fine.
        if (fine_exponent_offset < 0. || std::fabs(fine_exponent_offset - std::round(fine_exponent_offset)) > 1e-12) {
            return mio::failure(mio::StatusCode::InvalidValue,
                                "fine_exponent_offset has to be a non-negative integer, but is " +
                                    std::to_string(fine_exponent_offset) + ".");
        }

        ScalarType dt_ide     = pow(10, -ide_exponent);
        ScalarType div_dt_ide = pow(10, ide_exponent);
        std::cout << "Simulation with dt=" << dt_ide << std::endl;

        // Define the smaller step size and t_switch depending on the target step size dt_ide.
        ScalarType fine_exponent = ide_exponent + fine_exponent_offset;
        ScalarType t_switch      = t0_ide + num_steps_until_switch * dt_ide;

        // If the forward finite difference scheme is used, S' at the first finite_difference_order time points is
        // computed with values of S at later time points. Thus, the first simulation needs at least
        // 2 * finite_difference_order - 1 time steps.
        size_t num_fine_steps = num_steps_until_switch * (size_t)std::round(std::pow(10, fine_exponent_offset));
        if (num_steps_until_switch > 0 && forward_fd && num_fine_steps < 2 * finite_difference_order - 1) {
            return mio::failure(mio::StatusCode::InvalidValue,
                                fmt::format("The simulation until t_switch has only {} time steps, but at least {} "
                                            "are needed for the forward finite difference scheme.",
                                            num_fine_steps, 2 * finite_difference_order - 1));
        }

        // Values that are used to initialize the simulation with the target step size dt_ide. Without simulation with
        // step size dt_fine, only the values of the ODE solution at t0_ide are used. If flows_init is empty, the flows
        // on the initial time window are computed in the simulation.
        mio::TimeSeries<ScalarType> init_populations = init_values_t0;
        mio::TimeSeries<ScalarType> flows_init((size_t)mio::isir::InfectionTransition::Count);

        if (num_steps_until_switch > 0) {
            // Simulate until t_switch with the smaller step size and save the result as CSV files.
            std::string csv_dir  = save_dir.empty() ? "./" : save_dir;
            std::string csv_file = fmt::format(
                "{}/ide_result_until_tswitch_dt=1e-{:.0f}_dtfine=1e-{:.0f}_tswitch={:.10g}_gregoryorder={}.csv",
                csv_dir, ide_exponent, fine_exponent, t_switch, gregory_order);
            std::string csv_file_flows = fmt::format(
                "{}/ide_result_until_tswitch_dt=1e-{:.0f}_dtfine=1e-{:.0f}_tswitch={:.10g}_gregoryorder={}_flows.csv",
                csv_dir, ide_exponent, fine_exponent, t_switch, gregory_order);
            BOOST_OUTCOME_TRY(simulate_ide_until_switch(fine_exponent, t_switch, gregory_order, finite_difference_order,
                                                        init_values_t0, TimeInfected, cont_freq, csv_file,
                                                        csv_file_flows, kahan, forward_fd));

            // Read in the values obtained with the smaller step size and use the whole time interval
            // [t0_ide, t_switch] of the first simulation as initial time window.
            ScalarType init_window = num_steps_until_switch * dt_ide;
            BOOST_OUTCOME_TRY(auto&& compartments_until_switch,
                              read_timeseries_from_csv(csv_file, (size_t)mio::isir::InfectionState::Count));
            BOOST_OUTCOME_TRY(auto&& init_populations_coarse,
                              get_init_values(compartments_until_switch, init_window, t_switch, dt_ide, div_dt_ide));
            init_populations = init_populations_coarse;

            if (use_fine_flows) {
                BOOST_OUTCOME_TRY(
                    auto&& flows_until_switch,
                    read_timeseries_from_csv(csv_file_flows, (size_t)mio::isir::InfectionTransition::Count));
                BOOST_OUTCOME_TRY(auto&& flows_init_coarse,
                                  get_init_values(flows_until_switch, init_window, t_switch, dt_ide, div_dt_ide));
                flows_init = flows_init_coarse;
            }
        }

        // Simulate until tmax with the target step size.
        std::vector<mio::TimeSeries<ScalarType>> result =
            run_ide(std::move(init_populations), std::move(flows_init), dt_ide, div_dt_ide, gregory_order,
                    finite_difference_order, tmax, TimeInfected, cont_freq, kahan, forward_fd);

        if (!save_dir.empty()) {
            // Save compartments.
            mio::TimeSeries<ScalarType> compartments = result[0];
            mio::TimeSeries<ScalarType> flows        = result[1];

            auto csv_result = compartments.export_csv(fmt::format("{}/ide_result.csv", save_dir), {}, ',', 12);

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

    // The ODE model is initialized at t0_ode with the values from params. The IDE model is initialized at t0_ide with
    // the values of the ODE solution at t0_ide, so that both solutions can be compared at the same time points. Thus,
    // t0_ide has to be in [t0_ode, tmax] and has to be a time point of the saved ODE solution, i.e. a multiple of
    // 10^(-saving_exponent) after t0_ode.
    ScalarType t0_ode                   = 0.;
    ScalarType t0_ide                   = 50.;
    std::vector<ScalarType> tmax_values = {t0_ide + 100.};

    bool kahan      = false;
    bool forward_fd = true;

    std::vector<size_t> finite_difference_orders = {4};

    std::vector<size_t> gregory_orders = {1, 2, 3};

    // Before t_switch, the IDE model is simulated with the smaller step size dt_fine = dt * 10^(-fine_exponent_offset)
    // where dt is the target step size. The entries of fine_exponent_offsets correspond to the entries of
    // ide_exponents. The switching time is given by t_switch = t0_ide + num_steps_until_switch * dt with
    // num_steps_until_switch = max(gregory_order - 1, 2 * finite_difference_order - 1). Thus, there are at least
    // 2 * finite_difference_order - 1 time steps until t_switch which are needed for the forward finite difference
    // scheme, also if fine_exponent_offset = 0, i.e. dt_fine = dt. The simulation until t_switch does not use an
    // initial time window but is only initialized with the values of the ODE solution at t0_ide. Its result is saved
    // as CSV file and read in to initialize the simulation with step size dt, i.e. [t0_ide, t_switch] is the initial
    // time window of the simulation with step size dt. With fine_exponent_offset = 0, this yields the same result as a
    // simulation that directly uses the step size dt.
    // If use_fine_flows, the flows (i.e. -S') computed with step size dt_fine are also used on the initial time window
    // of the simulation with step size dt. Otherwise, they are computed from the compartments with step size dt.
    std::vector<ScalarType> ide_exponents         = {0., 1., 2., 3.};
    std::vector<ScalarType> fine_exponent_offsets = {0., 0., 0., 0.};
    bool use_fine_flows                           = true;

    // String with the entries of fine_exponent_offsets for the folder name, e.g. "1-1-2".
    std::string fine_exponent_offsets_string = "";
    for (size_t i = 0; i < fine_exponent_offsets.size(); i++) {
        fine_exponent_offsets_string += (i > 0 ? "-" : "") + fmt::format("{}", fine_exponent_offsets[i]);
    }

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
                fmt::format("./simulation_results/2026-10-16/"
                            "dtode=1e-{}_t0ode={}_timeinf={}_contfreq={}_kahan={}_forwardfd={}/"
                            "detailed_init_exponential_t0ide={}_tmax={}_finite_diff={}_fineoffsets={}_fineflows={}/",
                            ode_exponent, t0_ode, time_infected, cont_freq, kahan, forward_fd, t0_ide, tmax,
                            finite_difference_order, fine_exponent_offsets_string, use_fine_flows);

            // Make folder if not existent yet.
            std::filesystem::path dir(save_dir);
            std::filesystem::create_directories(dir);

            // ScalarType saving_exponent = *std::max_element(ide_exponents.begin(), ide_exponents.end());
            ScalarType saving_exponent = 3.;
            auto result_ode =
                simulate_ode(ode_exponent, t0_ode, tmax, time_infected, cont_freq, save_dir, saving_exponent).value();

            // The results of the simulation with step size dt start at t0_ide, hence we use t_init = t0_ide for the
            // folder name.
            std::string save_dir_ide = fmt::format("{}/tinit={}/", save_dir, t0_ide);
            // Make folder if not existent yet.
            std::filesystem::path dir_ide(save_dir_ide);
            std::filesystem::create_directories(dir_ide);

            // Do IDE simulations.
            for (size_t gregory_order : gregory_orders) {
                std::cout << std::endl;
                std::cout << "Gregory order: " << gregory_order << std::endl;
                // Number of time steps with step size dt that are simulated with the smaller step size dt_fine.
                size_t num_steps_until_switch = std::max(gregory_order - 1, 2 * finite_difference_order - 1);
                mio::IOResult<void> result_ide =
                    simulate_ide(ide_exponents, fine_exponent_offsets, gregory_order, finite_difference_order,
                                 result_ode[0], t0_ide, tmax, time_infected, cont_freq, num_steps_until_switch,
                                 use_fine_flows, save_dir_ide, kahan, forward_fd);
                if (!result_ide) {
                    std::cout << result_ide.error().formatted_message() << std::endl;
                    return -1;
                }
            }
        }
    }
}
