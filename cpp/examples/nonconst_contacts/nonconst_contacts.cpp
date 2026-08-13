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

ScalarType cont_freq = 0.73;

ScalarType S0               = 9999000.;
ScalarType I0               = 1000.;
ScalarType R0               = 0.;
ScalarType total_population = S0 + I0 + R0;
} // namespace params

mio::UncertainContactMatrix<ScalarType> scale_contact_matrix(ScalarType damping, ScalarType damping_time)
{
    using namespace params;

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup<ScalarType>(1, 1);

    // Perform simulation with a decrease in contacts.
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    // contact_matrix[0].add_damping(0., mio::SimulationTime<ScalarType>(0.));
    contact_matrix[0].add_damping(damping, mio::SimulationTime<ScalarType>(damping_time));

    return mio::UncertainContactMatrix(contact_matrix);
}

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

mio::IOResult<mio::TimeSeries<ScalarType>> simulate_ode(ScalarType ode_exponent, ScalarType t0_ode, ScalarType tmax,
                                                        ScalarType TimeInfected, ScalarType damping,
                                                        ScalarType damping_time, std::string save_dir = "",
                                                        ScalarType saving_exponent = 0.,
                                                        ScalarType smoother_window = 2., size_t smoothstep_order = 4)
{
    using namespace params;

    ScalarType dt_ode = pow(10, -ode_exponent);

    mio::log_info("Simulating ODE-SIR; t={} ... {} with dt = {}.", t0_ode, tmax, dt_ode);

    mio::osir::Model<ScalarType> model(num_agegroups, smoother_window, smoothstep_order);

    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Susceptible}] = S0;
    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Infected}]    = I0;
    model.populations[{mio::AgeGroup(0), mio::osir::InfectionState::Recovered}]   = R0;

    model.parameters.set<mio::osir::TimeInfected<ScalarType>>(TimeInfected);
    model.parameters.set<mio::osir::TransmissionProbabilityOnContact<ScalarType>>(TransmissionProbabilityOnContact);

    // mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    // contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));

    mio::UncertainContactMatrix<ScalarType> contact_matrix         = scale_contact_matrix(damping, damping_time);
    model.parameters.get<mio::osir::ContactPatterns<ScalarType>>() = mio::UncertainContactMatrix(contact_matrix);

    model.check_constraints();

    std::unique_ptr<mio::OdeIntegratorCore<ScalarType>> integrator =
        std::make_unique<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_fehlberg78>>();
    // std::unique_ptr<mio::OdeIntegratorCore<ScalarType>> integrator =
    //     std::make_unique<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_fehlberg78>>(
    //         1e-10, 1e-7, dt_ode, dt_ode);

    auto sir = simulate<ScalarType, mio::osir::Model<ScalarType>>(t0_ode, tmax, dt_ode, model, std::move(integrator));

    std::cout << "Num tps ODE: " << sir.get_num_time_points() << std::endl;

    mio::TimeSeries<ScalarType> compressed_compartments =
        mio::TimeSeries<ScalarType>((Eigen::Index)mio::isir::InfectionState::Count);

    if (!save_dir.empty()) {
        // Save compartments.
        compressed_compartments = compress_timeseries(sir, saving_exponent);
        std::cout << "Num tps ODE compressed: " << compressed_compartments.get_num_time_points() << std::endl;
        // auto result             = compartments.export_csv("ode_result.csv");
        auto save_result_status_ode =
            mio::save_result({compressed_compartments}, {0}, num_agegroups,
                             save_dir + "result_ode_dt=1e-" + fmt::format("{:.0f}", ode_exponent) + "_savedt=1e-" +
                                 fmt::format("{:.0f}", saving_exponent) + ".h5");

        if (!save_result_status_ode) {
            return mio::failure(mio::StatusCode::InvalidValue,
                                "Error occured while saving the ODE simulation results.");
        }
    }
    return mio::success(compressed_compartments);
}

mio::IOResult<void> simulate_ide(ScalarType ide_exponent, ScalarType ode_exponent, size_t gregory_order,
                                 size_t finite_difference_order, ScalarType t_init_window, ScalarType t0_ide,
                                 ScalarType tmax, ScalarType TimeInfected, ScalarType damping, ScalarType damping_time,
                                 ScalarType smoother_window = 2., size_t smoothstep_order = 4.,
                                 std::string save_dir = "",
                                 mio::TimeSeries<ScalarType> compartments_groundtruth =
                                     mio::TimeSeries<ScalarType>((size_t)mio::isir::InfectionState::Count),
                                 size_t fd_order_contacts = 4)
{
    using namespace params;
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType dt_ide     = pow(10, -ide_exponent);
    ScalarType div_dt_ide = pow(10, ide_exponent);

    std::cout << "Simulation with dt=" << dt_ide << std::endl;

    mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);

    if (compartments_groundtruth.get_num_time_points() == 0) {
        std::cout << "No groundtruth was given.\n";
    }
    else {
        std::cout << "Initializing with given groundtruth.\n";

        // Initialize time points before t0_ide based on groundtruth.
        ScalarType div_dt_groundtruth       = std::pow(10, ode_exponent);
        ScalarType groundtruth_index_factor = std::pow(10, ode_exponent - ide_exponent);

        Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));

        std::vector<size_t> compartments = {(size_t)mio::isir::InfectionState::Susceptible,
                                            (size_t)mio::isir::InfectionState::Infected,
                                            (size_t)mio::isir::InfectionState::Recovered};

        ScalarType t_init = t0_ide - t_init_window;
        int start_index   = std::round(t_init * div_dt_groundtruth);

        // Add values to init_populations.
        for (size_t compartment : compartments) {
            vec_init[compartment] = compartments_groundtruth.get_value(start_index)[compartment];
        }

        init_populations.add_time_point(t_init, vec_init);

        while (init_populations.get_last_time() < t0_ide - 1e-10) {
            std::round(start_index += groundtruth_index_factor);
            for (size_t compartment : compartments) {
                vec_init[compartment] = compartments_groundtruth.get_value(std::round(start_index))[compartment];
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
    // mio::InitialZeroInfectiousness transmissiononcontact(1., 5., TransmissionProbabilityOnContact);
    mio::StateAgeFunctionWrapper transmissiononcontact_wrapper(transmissiononcontact);
    model.parameters.get<mio::isir::TransmissionProbabilityOnContact>() = transmissiononcontact_wrapper;

    mio::ConstantFunction riskofinfection(RiskOfInfectionFromSymptomatic);
    mio::StateAgeFunctionWrapper riskofinfection_wrapper(riskofinfection);
    model.parameters.get<mio::isir::RiskOfInfectionFromSymptomatic>() = riskofinfection_wrapper;

    // mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup(1, 1);
    // contact_matrix[0]                      = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    mio::UncertainContactMatrix<ScalarType> contact_matrix = scale_contact_matrix(damping, damping_time);
    model.parameters.get<mio::isir::ContactPatterns>()     = mio::UncertainContactMatrix(contact_matrix);

    std::cout << "support max: " << model.compute_calctime(dt_ide, 1e-6) << std::endl;

    // Carry out simulation.
    mio::isir::SimulationMessinaExtendedDetailedInit sim(model, dt_ide, div_dt_ide);
    sim.advance(tmax, true, fd_order_contacts, damping_time, smoother_window, smoothstep_order);

    if (!save_dir.empty()) {
        // Save compartments.
        mio::TimeSeries<ScalarType> compartments = sim.get_result();
        auto save_result_status_ide =
            mio::save_result({compartments}, {0}, num_agegroups,
                             save_dir + "result_ide_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
                                 +"_gregoryorder=" + fmt::format("{}", gregory_order) + ".h5");

        if (!save_result_status_ide) {
            return mio::failure(mio::StatusCode::InvalidValue,
                                "Error occured while saving the IDE simulation results.");
        }
    }

    return mio::success();
}

int main()
{
    using namespace params;

    ScalarType time_infected = 2.;

    ScalarType t0_ode      = 0.;
    ScalarType t0_ide      = 50.;
    ScalarType init_window = 40.;
    ScalarType tmax        = 60.;

    ScalarType damping      = 0.1;
    ScalarType damping_time = 55.372907;

    std::vector<size_t> gregory_orders = {1, 2, 3};
    size_t finite_difference_order     = 4;
    size_t fd_order_contacts           = 4;
    size_t smoothstep_order            = 0; // possible values: 0, 3 or 4; smoothstep_order=0 leads to smoother_cosine
    ScalarType smoother_window         = 2.;

    // Compute groundtruth with ODE model.
    ScalarType ode_exponent               = 6.;
    std::vector<ScalarType> ide_exponents = {0, 1, 2};

    std::string save_dir =
        fmt::format("./simulation_results/2026-08-13/"
                    "phi_deriv_analytical_smoothstepc{}_fdordercontacts={}_smootherwindow={}_t0ode={}/"
                    "nonconst_contacts_t0ide={}_tmax={}_dampingtime={}_damping={}/",
                    smoothstep_order, fd_order_contacts, smoother_window, t0_ode, t0_ide, tmax, damping_time, damping);

    // Make folder if not existent yet.
    std::filesystem::path dir(save_dir);
    std::filesystem::create_directories(dir);

    // ScalarType saving_exponent = *std::max_element(ide_exponents.begin(), ide_exponents.end());
    ScalarType saving_exponent = 3.;
    std::cout << "saving exp: " << saving_exponent << std::endl;
    auto result_ode = simulate_ode(ode_exponent, t0_ode, tmax, time_infected, damping, damping_time, save_dir,
                                   saving_exponent, smoother_window, smoothstep_order)
                          .value();

    ScalarType t_init        = t0_ide - init_window;
    std::string save_dir_ide = fmt::format("{}/tinit={}/", save_dir, t_init);
    // Make folder if not existent yet.
    std::filesystem::path dir_ide(save_dir_ide);
    std::filesystem::create_directories(dir_ide);

    // Do IDE simulations.
    for (size_t gregory_order : gregory_orders) {
        std::cout << "Gregory order: " << gregory_order << std::endl;
        for (ScalarType ide_exponent : ide_exponents) {
            std::cout << std::endl;
            mio::IOResult<void> result_ide =
                simulate_ide(ide_exponent, saving_exponent, gregory_order, finite_difference_order, init_window, t0_ide,
                             tmax, time_infected, damping, damping_time, smoother_window, smoothstep_order,
                             save_dir_ide, result_ode, fd_order_contacts);
        }
    }
}