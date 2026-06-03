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

ScalarType TransmissionProbabilityOnContact = 0.5;
ScalarType RiskOfInfectionFromSymptomatic   = 1.;
ScalarType Seasonality                      = 0.;

ScalarType cont_freq = 1.;

ScalarType S0               = 999000.;
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

mio::IOResult<mio::TimeSeries<ScalarType>>
simulate_ide(ScalarType ide_exponent, size_t gregory_order, size_t finite_difference_order, ScalarType t_init,
             ScalarType t0, ScalarType tmax, ScalarType TimeInfected, ScalarType damping, ScalarType damping_time,
             std::string save_dir = "", std::string model_type = "",
             mio::TimeSeries<ScalarType> result_groundtruth =
                 mio::TimeSeries<ScalarType>((size_t)mio::isir::InfectionState::Count))
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
        //         S0 - 500 * init_populations.get_num_time_points() * dt * dt;
        //     vec_init[(size_t)mio::isir::InfectionState::Infected] =
        //         I0 + 500 * init_populations.get_num_time_points() * dt * dt;
        //     vec_init[(size_t)mio::isir::InfectionState::Recovered] = R0;

        //     init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
        // }
    }
    else {
        std::cout << "Initializing with given groundtruth for t_{-4},...,t0.\n";

        if (model_type == "simple") {

            // Initialize time points before t0 based on groundtruth.
            // Get index of t0 in groundtruth.
            size_t t0_index = 0;
            for (Eigen::Index i = 0; i < result_groundtruth.get_num_time_points(); i++) {
                ScalarType t = result_groundtruth.get_time(i);
                if (fabs(t - t0) < 1e-7) {
                    t0_index = i;
                }
            }
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

    mio::ExponentialSurvivalFunction logn(1. / TimeInfected);

    mio::StateAgeFunctionWrapper dist(logn);
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

    std::cout << "support max: " << model.compute_calctime(dt, 1e-6) << std::endl;

    // Carry out simulation.
    mio::isir::SimulationMessinaExtendedDetailedInit sim(model, dt);
    sim.advance(tmax);
    mio::TimeSeries<ScalarType> compartments = sim.get_result();

    if (!save_dir.empty()) {
        // Save compartments.

        auto save_result_status_ide =
            mio::save_result({compartments}, {0}, num_agegroups,
                             save_dir + model_type + "_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
                                 +"_gregoryorder=" + fmt::format("{}", gregory_order) + ".h5");
    }

    return mio::success(compartments);
}

int main()
{
    using namespace params;

    ScalarType time_infected = 5.;

    ScalarType t_init = 0.;

    ScalarType t0              = 20.;
    ScalarType t_init_simple   = t0;
    ScalarType t_init_detailed = t0 - 5.;

    ScalarType tmax = t0 + 20.;

    ScalarType damping      = 0.;
    ScalarType damping_time = 0.;

    size_t gregory_order           = 3;
    size_t finite_difference_order = 4;

    ScalarType ide_exponent = 2.;

    std::string save_dir = fmt::format("../../simulation_results/2026-06-01/diff_init_conditions_timeinf={}/"
                                       "nonconst_contacts_tinit={}_t0={}_tmax={}/",
                                       time_infected, t_init, t0, tmax);

    // Make folder if not existent yet.
    std::filesystem::path dir(save_dir);
    std::filesystem::create_directories(dir);

    // Do IDE simulations.
    std::cout << std::endl;
    mio::IOResult<mio::TimeSeries<ScalarType>> result_ide_groundtruth =
        simulate_ide(ide_exponent, gregory_order, finite_difference_order, t_init, t_init, tmax, time_infected, damping,
                     damping_time, save_dir, "groundtruth");

    mio::IOResult<mio::TimeSeries<ScalarType>> result_ide_detailed =
        simulate_ide(ide_exponent, gregory_order, finite_difference_order, t_init_detailed, t0, tmax, time_infected,
                     damping, damping_time, save_dir, "detailed", result_ide_groundtruth.value());

    mio::IOResult<mio::TimeSeries<ScalarType>> result_ide_simple =
        simulate_ide(ide_exponent, gregory_order, finite_difference_order, t_init_simple, t0, tmax, time_infected,
                     damping, damping_time, save_dir, "simple", result_ide_groundtruth.value());
}
