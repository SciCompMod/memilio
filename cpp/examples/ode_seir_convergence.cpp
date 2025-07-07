/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/type_list.h"
#include "state_estimators.h"
#include <cmath> // Required for std::ceil

using namespace mio::examples;
using SimType = mio::FlowSimulation<ScalarType, mio::oseir::Model<ScalarType>>;

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0;
    ScalarType tmax = 1.;

    mio::oseir::Model<ScalarType> model(1);

    const auto sus = 9700, exp = 100, inf = 100, rec = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = sus;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = exp;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]    = inf;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}]   = rec;

    model.parameters.set<mio::oseir::TimeExposed<ScalarType>>(5.2);
    model.parameters.set<mio::oseir::TimeInfected<ScalarType>>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.1);

    // mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();
    // contact_matrix[0].get_baseline().setConstant(2.7);
    // contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.));

    model.check_constraints();

    // Total solution
    auto integrator = std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
    const double step_size_ref = 1e-8;
    auto sim                   = mio::FlowSimulation<double, mio::oseir::Model<double>>(model, t0, step_size_ref);
    sim.set_integrator(integrator);
    sim.advance(tmax);
    const auto& seir_res = sim.get_result();

    // commuter state estimations
    const auto p_c   = 0.3;
    const double S_c = sus * p_c;
    const double E_c = exp * p_c;
    const double I_c = inf * p_c;
    const double R_c = rec * p_c;
    Eigen::VectorXd initial_mobile_pop(4);
    initial_mobile_pop << S_c, E_c, I_c, R_c;

    // calculate high-precision reference solution for mobile population using ExplicitModel
    mio::log_info("Calculating high-precision reference sol using ExplicitModel...");
    ExplicitModel model_explicit_ref(1, 1); // 1 age group, 1 commuter group (CommuterBase)
    const ScalarType p_nc = 1.0 - p_c;

    // Initialize populations for the explicit reference model
    // NonCommuter
    model_explicit_ref.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::S}] = sus * p_nc;
    model_explicit_ref.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::E}] = exp * p_nc;
    model_explicit_ref.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::I}] = inf * p_nc;
    model_explicit_ref.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::NonCommuter, InfectionStateExplicit::R}] = rec * p_nc;

    // Commuter (using CommuterBase for the single commuter group)
    model_explicit_ref.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::S}] = S_c;
    model_explicit_ref.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::E}] = E_c;
    model_explicit_ref.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::I}] = I_c;
    model_explicit_ref.populations[mio::Index<mio::AgeGroup, CommuterType, InfectionStateExplicit>{
        mio::AgeGroup(0), CommuterType::CommuterBase, InfectionStateExplicit::R}] = R_c;

    // Copy parameters from the aggregate flow model to the explicit reference model
    model_explicit_ref.parameters.set<mio::oseir::TimeExposed<ScalarType>>(
        model.parameters.get<mio::oseir::TimeExposed<ScalarType>>()[mio::AgeGroup(0)]);
    model_explicit_ref.parameters.set<mio::oseir::TimeInfected<ScalarType>>(
        model.parameters.get<mio::oseir::TimeInfected<ScalarType>>()[mio::AgeGroup(0)]);
    model_explicit_ref.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(
        model.parameters.get<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>()[mio::AgeGroup(0)]);
    model_explicit_ref.parameters.get<mio::oseir::ContactPatterns<ScalarType>>() =
        model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();

    model_explicit_ref.check_constraints();

    auto integrator_ref =
        std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
    ExplicitSim sim_explicit_ref(model_explicit_ref, t0, step_size_ref);
    sim_explicit_ref.set_integrator(integrator_ref);
    sim_explicit_ref.advance(tmax);
    const auto& explicit_ref_results_full = sim_explicit_ref.get_result();

    // The commuter population is the last 4 states in the explicit model with 1 age group and 1 commuter type
    auto y_ref_end = explicit_ref_results_full.get_last_value().tail(4);

    // Convergence test setup
    const std::vector<ScalarType> dt_values = {1e-0, 0.5, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
    struct ConvergenceResult {
        ScalarType dt;
        int time_steps; // Added to store the actual number of steps
        ScalarType error_euler;
        ScalarType error_high_order;
        ScalarType error_flow;
        ScalarType error_prob;
    };
    std::vector<ConvergenceResult> convergence_data;

    mio::log_info("Starting convergence test!");

    for (ScalarType current_dt : dt_values) {
        mio::log_info("dt = {}", current_dt);
        Eigen::VectorXd mobile_pop_temp;

        // Use std::ceil to ensure t0 + num_steps * current_dt >= tmax.
        const int num_steps = static_cast<int>(std::ceil((tmax - t0) / current_dt));
        // const int num_steps = static_cast<int>((tmax - t0) / current_dt + 1e-9); // Alternative way to calculate num_steps

        // Euler
        mobile_pop_temp      = initial_mobile_pop;
        ScalarType t_current = t0;
        for (int step = 0; step < num_steps; ++step) {
            const auto closest_idx_total = static_cast<size_t>(std::round(t_current / step_size_ref));
            const auto& total_pop_at_t   = seir_res.get_value(closest_idx_total);
            integrate_mobile_population_euler(mobile_pop_temp, sim, total_pop_at_t, t_current, current_dt);
            t_current += current_dt;
        }
        const auto y_euler_end       = mobile_pop_temp;
        const ScalarType error_euler = (y_euler_end - y_ref_end).norm();

        // high Order
        mobile_pop_temp = initial_mobile_pop;
        t_current       = t0;
        for (int step = 0; step < num_steps; ++step) {
            const auto closest_idx_total = static_cast<size_t>(std::round(t_current / step_size_ref));
            const auto& total_pop_at_t   = seir_res.get_value(closest_idx_total);
            integrate_mobile_population_high_order(mobile_pop_temp, sim, total_pop_at_t, t_current, current_dt);
            t_current += current_dt;
        }
        const auto y_high_order_end       = mobile_pop_temp;
        const ScalarType error_high_order = (y_high_order_end - y_ref_end).norm();

        // Flow-based
        mobile_pop_temp = initial_mobile_pop;
        t_current       = t0;
        for (int step = 0; step < num_steps; ++step) {
            const auto closest_idx_total = static_cast<size_t>(std::round(t_current / step_size_ref));
            const auto& total_pop_at_t   = seir_res.get_value(closest_idx_total);
            flow_based_mobility_returns(mobile_pop_temp, sim, total_pop_at_t, t_current, current_dt);
            t_current += current_dt;
        }
        const auto y_flow_end       = mobile_pop_temp;
        const ScalarType error_flow = (y_flow_end - y_ref_end).norm();

        // Probabilistic
        mobile_pop_temp = initial_mobile_pop;
        t_current       = t0;
        for (int step = 0; step < num_steps; ++step) {
            const auto closest_idx_total = static_cast<size_t>(std::round(t_current / step_size_ref));
            const auto& total_pop_at_t   = seir_res.get_value(closest_idx_total);
            probabilistic_mobility_returns(mobile_pop_temp, sim, total_pop_at_t, t_current, current_dt);
            t_current += current_dt;
        }
        const auto y_prob_end       = mobile_pop_temp;
        const ScalarType error_prob = (y_prob_end - y_ref_end).norm();

        // Store results
        convergence_data.push_back({current_dt, num_steps, error_euler, error_high_order, error_flow, error_prob});
    }
    // Write to file
    const std::string save_dir             = "/localdata1/code/memilio/saves";
    const std::string convergence_filename = save_dir + "/convergence_results.csv";
    std::ofstream convergence_file(convergence_filename);
    convergence_file << "dt,TimeSteps,ErrorEuler,ErrorHighOrder,ErrorFlow,ErrorProb\n";
    for (const auto& result : convergence_data) {
        convergence_file << result.dt << "," << result.time_steps << "," << result.error_euler << ","
                         << result.error_high_order << "," << result.error_flow << "," << result.error_prob << "\n";
    }
    convergence_file.close();

    return 0;
}
