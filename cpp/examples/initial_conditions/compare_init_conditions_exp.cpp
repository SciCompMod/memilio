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
#include "memilio/compartments/flow_simulation.h"
#include "memilio/epidemiology/contact_matrix.h"
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

ScalarType TransmissionProbabilityOnContact = 0.8;
ScalarType RiskOfInfectionFromSymptomatic   = 1.;
ScalarType Seasonality                      = 0.;

ScalarType TimeInfected = 2.;

ScalarType cont_freq = 0.73;

ScalarType total_population = 1e7;
ScalarType I0               = 1000.;
ScalarType R0               = 0.;
ScalarType S0               = total_population - I0 - R0;
} // namespace params

mio::IOResult<mio::TimeSeries<ScalarType>> simulate_ode(ScalarType ode_exponent, ScalarType t0_ode, ScalarType tmax,
                                                        std::string save_dir = "")
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
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, cont_freq));
    model.parameters.get<mio::osir::ContactPatterns<ScalarType>>() = mio::UncertainContactMatrix(contact_matrix);

    model.check_constraints();

    std::unique_ptr<mio::OdeIntegratorCore<ScalarType>> integrator =
        std::make_unique<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_fehlberg78>>();

    auto sim = mio::FlowSimulation<ScalarType, mio::osir::Model<ScalarType>>(model, t0_ode, dt_ode);
    sim.set_integrator_core(std::move(integrator));
    sim.set_last_step_tolerance(1e-8);

    sim.advance(tmax);
    auto compartments = sim.get_result();

    std::cout << "Num tps ODE: " << compartments.get_num_time_points() << std::endl;

    if (!save_dir.empty()) {
        // Save compartments.

        auto result = compartments.export_csv(fmt::format("{}/ode_result_compressed.csv", save_dir));
        auto save_result_status_ode =
            mio::save_result({compartments}, {0}, num_agegroups,
                             save_dir + "result_ode_dt=1e-" + fmt::format("{:.0f}", ode_exponent) + "_savedt=1e-" +
                                 fmt::format("{:.0f}", ode_exponent) + ".h5");

        if (!save_result_status_ode) {
            return mio::failure(mio::StatusCode::InvalidValue,
                                "Error occured while saving the ODE simulation results.");
        }
    }

    return mio::success(compartments);
}

mio::IOResult<mio::TimeSeries<ScalarType>>
simulate_ide(ScalarType ide_exponent, size_t gregory_order, size_t finite_difference_order, ScalarType t_init,
             ScalarType t0, ScalarType tmax, std::string save_dir = "", std::string model_type = "",
             mio::TimeSeries<ScalarType> result_groundtruth =
                 mio::TimeSeries<ScalarType>((size_t)mio::isir::InfectionState::Count),
             size_t write_inf_per_infage_time = 0)
{
    using namespace params;
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType dt = pow(10, -ide_exponent);
    std::cout << "Simulation with dt=" << dt << std::endl;

    mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);

    if (result_groundtruth.get_num_time_points() == 0) {
        std::cout << "No groundtruth was given.\n";
        Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));

        vec_init[(size_t)mio::isir::InfectionState::Susceptible] = S0;
        vec_init[(size_t)mio::isir::InfectionState::Infected]    = I0;
        vec_init[(size_t)mio::isir::InfectionState::Recovered]   = R0;

        init_populations.add_time_point(t_init, vec_init);

        // while (init_populations.get_last_time() < t0 - 1e-10) {

        //     vec_init[(size_t)mio::isir::InfectionState::Susceptible] =
        //         S0 - 5000 * init_populations.get_num_time_points() * dt;
        //     vec_init[(size_t)mio::isir::InfectionState::Infected] =
        //         I0 + 5000 * init_populations.get_num_time_points() * dt;
        //     vec_init[(size_t)mio::isir::InfectionState::Recovered] =
        //         R0 + 4000 * init_populations.get_num_time_points() * dt;

        //     init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
        // }
    }
    else {
        if (model_type == "simple") {

            // Initialize time points before t0 based on groundtruth.
            // Get index of t0 in groundtruth.
            size_t t0_index = size_t((t0 - result_groundtruth.get_time(0)) / dt);
            std::cout << "t0 index: " << t0_index << std::endl;

            Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));

            vec_init[(size_t)mio::isir::InfectionState::Susceptible] =
                result_groundtruth.get_value(t0_index)[(size_t)mio::isir::InfectionState::Susceptible];
            vec_init[(size_t)mio::isir::InfectionState::Infected] =
                result_groundtruth.get_value(t0_index)[(size_t)mio::isir::InfectionState::Infected];
            vec_init[(size_t)mio::isir::InfectionState::Recovered] =
                result_groundtruth.get_value(t0_index)[(size_t)mio::isir::InfectionState::Recovered];

            init_populations.add_time_point(t0, vec_init);
        }

        else {

            Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));

            std::vector<size_t> compartments = {(size_t)mio::isir::InfectionState::Susceptible,
                                                (size_t)mio::isir::InfectionState::Infected,
                                                (size_t)mio::isir::InfectionState::Recovered};

            // Add values to init_populations.
            for (size_t compartment : compartments) {
                vec_init[compartment] = result_groundtruth.get_value(t_init / dt)[compartment];
            }
            std::cout << "groundtruth num tps: " << result_groundtruth.get_num_time_points() << std::endl;
            std::cout << t_init / dt << std::endl;
            std::cout << "t_init: " << result_groundtruth.get_time(t_init / dt) << std::endl;

            init_populations.add_time_point(t_init, vec_init);

            while (init_populations.get_last_time() < t0 - 1e-10) {
                for (size_t compartment : compartments) {
                    vec_init[compartment] =
                        result_groundtruth.get_value(t_init / dt + init_populations.get_num_time_points())[compartment];
                }
                init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
            }
        }

        // std::cout << "Init populations tps:  \n";
        // for (Eigen::Index i = 0; i < init_populations.get_num_time_points(); i++) {
        //     std::cout << init_populations.get_time(i) << std::endl;
        // }
    }

    // Initialize model.
    mio::isir::ModelMessinaExtendedDetailedInit model(std::move(init_populations), total_population, gregory_order,
                                                      finite_difference_order);

    mio::ExponentialSurvivalFunction survival_func(1. / TimeInfected);
    std::cout << "mean: " << survival_func.get_mean(dt) << std::endl;
    mio::StateAgeFunctionWrapper dist(survival_func);

    std::vector<mio::StateAgeFunctionWrapper<ScalarType>> vec_dist((size_t)mio::isir::InfectionTransition::Count, dist);
    model.parameters.get<mio::isir::TransitionDistributions>() = vec_dist;

    mio::ConstantFunction transmissiononcontact(TransmissionProbabilityOnContact);
    // mio::InitialZeroInfectiousness transmissiononcontact(1., 5., TransmissionProbabilityOnContact);
    mio::StateAgeFunctionWrapper transmissiononcontact_wrapper(transmissiononcontact);
    model.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = transmissiononcontact_wrapper;

    mio::ConstantFunction riskofinfection(RiskOfInfectionFromSymptomatic);
    mio::StateAgeFunctionWrapper riskofinfection_wrapper(riskofinfection);
    model.parameters.get<mio::isir::RiskOfInfectionFromSymptomatic>() = riskofinfection_wrapper;

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup<ScalarType>(1, 1);
    contact_matrix[0]                      = mio::ContactMatrix<ScalarType>(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    model.parameters.get<mio::isir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    std::cout << "support max: " << model.compute_calctime(dt, 1e-6) << std::endl;

    // Carry out simulation.
    mio::isir::SimulationMessinaExtendedDetailedInit sim(model, dt);
    sim.advance(tmax);
    mio::TimeSeries<ScalarType> compartments = sim.get_result();
    mio::TimeSeries<ScalarType> flows        = sim.get_flows();

    if (!save_dir.empty()) {
        // Save compartments.

        auto save_result_status_ide =
            mio::save_result({compartments}, {0}, num_agegroups,
                             save_dir + model_type + "_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
                                 +"_gregoryorder=" + fmt::format("{}", gregory_order) + ".h5");

        auto save_result_status_ide_flows =
            mio::save_result({flows}, {0}, num_agegroups,
                             save_dir + model_type + "_flows_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
                                 +"_gregoryorder=" + fmt::format("{}", gregory_order) + ".h5");

        std::cout << "write time: " << write_inf_per_infage_time << std::endl;
        if (write_inf_per_infage_time > 0) {

            size_t write_inf_per_infage_index = (write_inf_per_infage_time - compartments.get_time(0)) / dt;
            std::cout << "write index: " << write_inf_per_infage_index << std::endl;
            std::cout << "eval time: " << compartments.get_time(write_inf_per_infage_index) << std::endl;
            std::cout << std::endl;
            mio::TimeSeries<ScalarType> infected_per_infection_age =
                sim.write_infected_per_infection_age(write_inf_per_infage_index);

            auto save_result_status_ide_inf = mio::save_result(
                {infected_per_infection_age}, {0}, num_agegroups,
                save_dir + model_type + "_infectionagedistribution_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
                    +"_gregoryorder=" + fmt::format("{}", gregory_order) + ".h5");
        }
    }

    return mio::success(compartments);
}

int main()
{
    using namespace params;

    // ScalarType time_infected = 2.;

    ScalarType t0_ode = 0.;
    ScalarType t0_ide = 50.;

    ScalarType t_init_simple         = t0_ide;
    ScalarType t_init_detailed       = t0_ide - 40.;
    ScalarType t_init_detailed_short = t0_ide - 10.;

    ScalarType tmax = 150.;

    size_t gregory_order           = 3;
    size_t finite_difference_order = 4;

    ScalarType dt_exponent = 2.; // Used for both ODE and IDE simulations

    std::string save_dir = fmt::format("../../simulation_results/2026-07-29/compare_different_inits_exp/"
                                       "nonconst_contacts_tinitgroundtruth={}_tinitdetailed={}_t0ide={}_tmax={}/",
                                       t0_ode, t_init_detailed, t0_ide, tmax);

    // Make folder if not existent yet.
    std::filesystem::path dir(save_dir);
    std::filesystem::create_directories(dir);

    // Do IDE simulations.
    std::cout << std::endl;
    mio::IOResult<mio::TimeSeries<ScalarType>> result_groundtruth = simulate_ode(dt_exponent, t0_ode, tmax, save_dir);

    mio::IOResult<mio::TimeSeries<ScalarType>> result_ide_detailed =
        simulate_ide(dt_exponent, gregory_order, finite_difference_order, t_init_detailed, t0_ide, tmax, save_dir,
                     "detailed", result_groundtruth.value(), t0_ide);

    mio::IOResult<mio::TimeSeries<ScalarType>> result_ide_detailed_short =
        simulate_ide(dt_exponent, gregory_order, finite_difference_order, t_init_detailed_short, t0_ide, tmax, save_dir,
                     "detailed_short", result_groundtruth.value(), t0_ide);

    mio::IOResult<mio::TimeSeries<ScalarType>> result_ide_simple =
        simulate_ide(dt_exponent, gregory_order, finite_difference_order, t_init_simple, t0_ide, tmax, save_dir,
                     "simple", result_groundtruth.value(), t0_ide);
}
