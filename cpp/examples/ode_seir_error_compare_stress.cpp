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
#include "state_estimators.h"

using namespace mio::examples;
using SimType = mio::FlowSimulation<ScalarType, mio::oseir::Model<ScalarType>>;

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0;
    ScalarType tmax = 0.5;

    mio::oseir::Model<ScalarType> model(1);

    const auto sus = 5000, exp = 1500, inf = 1500, rec = 2000;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = sus;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = exp;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]    = inf;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}]   = rec;

    model.parameters.set<mio::oseir::TimeExposed<ScalarType>>(1.);
    model.parameters.set<mio::oseir::TimeInfected<ScalarType>>(1.);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.1);

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::oseir::ContactPatterns<ScalarType>>();
    contact_matrix[0].get_baseline().setConstant(2.7);
    // contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.));

    model.check_constraints();

    // Total solution
    auto integrator = std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();

    const double step_size_ref = 1e-1;
    auto sim                   = mio::FlowSimulation<double, mio::oseir::Model<double>>(model, t0, step_size_ref);
    sim.set_integrator(integrator);
    sim.advance(tmax);
    const auto& seir_res = sim.get_result();

    // auto seir = mio::simulate_flows<double, mio::oseir::Model<double>>(t0, tmax, dt, model, integrator);

    // commuter state estimations
    const auto p_c   = 0.98;
    const double S_c = sus * p_c;
    const double E_c = exp * p_c;
    const double I_c = inf * p_c;
    const double R_c = rec * p_c;
    Eigen::VectorXd mobile_population(4);
    mobile_population << S_c, E_c, I_c, R_c;

    // calculate reference solution
    const auto sol_ref = calculate_mobile_population(t0, tmax, step_size_ref, flow_based_mobility_returns, sim,
                                                     seir_res, mobile_population);

    // Solutions for Euler
    const ScalarType common_step_size = 0.5;
    mobile_population << S_c, E_c, I_c, R_c;
    const auto sol_euler = calculate_mobile_population(t0, tmax, common_step_size, integrate_mobile_population_euler,
                                                       sim, seir_res, mobile_population);

    // Write to csv files
    const std::string save_dir                        = "/localdata1/code/memilio/saves";
    const std::vector<std::string> compartment_labels = {"S", "E", "I", "R"};
    sol_ref.export_csv(save_dir + "/mobile_ref_solution_stress.csv", compartment_labels);
    sol_euler.export_csv(save_dir + "/mobile_euler_solution_stress.csv", compartment_labels);
    seir_res.export_csv(save_dir + "/seir_solution_stress.csv", compartment_labels);

    return 0;
}
