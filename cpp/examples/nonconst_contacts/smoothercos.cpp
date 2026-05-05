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

#include "ide_sir/model_smoothercos.h"
#include "ide_sir/infection_state.h"
#include "ide_sir/parameters.h"
#include "ide_sir/simulation_smoothercos.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/utils/time_series.h"
#include "memilio/io/result_io.h"
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <dirent.h>
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

ScalarType smoothercos_deriv(ScalarType current_time, mio::UncertainContactMatrix<ScalarType> contact_matrix,
                             ScalarType damping_time, ScalarType smoother_window)
{
    ScalarType xleft  = damping_time - smoother_window;
    ScalarType xright = damping_time;

    if (current_time <= xleft || current_time >= xright) {
        return 0.;
    }

    ScalarType yleft  = contact_matrix.get_cont_freq_mat().get_matrix_at(SimulationTime<ScalarType>(xleft))(0, 0);
    ScalarType yright = contact_matrix.get_cont_freq_mat().get_matrix_at(SimulationTime<ScalarType>(xright))(0, 0);

    ScalarType deriv = -0.5 * (yleft - yright) / (xright - xleft) * std::numbers::pi_v<ScalarType> *
                       std::sin(std::numbers::pi_v<ScalarType> / (xright - xleft) * (current_time - xleft));

    return deriv;
}

ScalarType smoothstep(ScalarType current_time, ScalarType damping_time, ScalarType damping, ScalarType cont_freq,
                      ScalarType smoother_window)
{
    ScalarType xleft  = damping_time - smoother_window;
    ScalarType xright = damping_time;

    ScalarType yleft  = cont_freq;
    ScalarType yright = (1. - damping) * cont_freq;

    if (current_time <= xleft) {
        return yleft;
    }
    if (current_time >= xright) {
        return yright;
    }

    else {
        ScalarType normalized_time = (current_time - xleft) / (xright - xleft);

        ScalarType smoothed_value =
            yleft + (yright - yleft) * (3. * std::pow(normalized_time, 2) - 2 * std::pow(normalized_time, 3));

        return smoothed_value;
    }
}

ScalarType smoothstep_deriv(ScalarType current_time, ScalarType damping_time, ScalarType damping, ScalarType cont_freq,
                            ScalarType smoother_window)
{
    ScalarType xleft  = damping_time - smoother_window;
    ScalarType xright = damping_time;

    ScalarType yleft  = cont_freq;
    ScalarType yright = (1. - damping) * cont_freq;

    ScalarType deriv = 0.;

    if (current_time <= xleft || current_time >= xright) {
        deriv = 0.;
    }
    else {
        ScalarType normalized_time       = (current_time - xleft) / (xright - xleft);
        ScalarType normalized_time_deriv = 1. / (xright - xleft);
        deriv                            = (yright - yleft) * (6. * normalized_time * normalized_time_deriv -
                                    6 * std::pow(normalized_time, 2) * normalized_time_deriv);
    }
    return deriv;
}

ScalarType smoothstep_c2(ScalarType current_time, ScalarType damping_time, ScalarType damping, ScalarType cont_freq,
                         ScalarType smoother_window)
{

    ScalarType xleft  = damping_time - smoother_window;
    ScalarType xright = damping_time;

    ScalarType yleft  = cont_freq;
    ScalarType yright = (1. - damping) * cont_freq;

    if (current_time <= xleft) {
        return yleft;
    }
    if (current_time >= xright) {
        return yright;
    }

    else {
        ScalarType normalized_time = (current_time - xleft) / (xright - xleft);

        ScalarType smoothed_value =
            yleft + (yright - yleft) * (6. * std::pow(normalized_time, 5) - 15. * std::pow(normalized_time, 4) +
                                        10. * std::pow(normalized_time, 3));

        return smoothed_value;
    }
}

ScalarType smoothstep_c2_deriv(ScalarType current_time, ScalarType damping_time, ScalarType damping,
                               ScalarType cont_freq, ScalarType smoother_window)
{
    ScalarType xleft  = damping_time - smoother_window;
    ScalarType xright = damping_time;

    ScalarType yleft  = cont_freq;
    ScalarType yright = (1. - damping) * cont_freq;

    ScalarType deriv = 0.;

    if (current_time <= xleft || current_time >= xright) {
        deriv = 0.;
    }
    else {
        ScalarType normalized_time       = (current_time - xleft) / (xright - xleft);
        ScalarType normalized_time_deriv = 1. / (xright - xleft);
        deriv                            = (yright - yleft) *
                (30. * std::pow(normalized_time, 4) - 60. * std::pow(normalized_time, 3) +
                 30 * std::pow(normalized_time, 2)) *
                normalized_time_deriv;
    }

    return deriv;
}

ScalarType sigmoid(ScalarType x)
{
    return 1. / (1. + std::exp(-x));
}

ScalarType sigmoid_smoother(ScalarType current_time, ScalarType damping_time, ScalarType damping, ScalarType cont_freq,
                            ScalarType smoother_window, ScalarType sigmoid_param)
{

    ScalarType xleft  = damping_time - smoother_window;
    ScalarType xright = damping_time;

    ScalarType yleft  = cont_freq;
    ScalarType yright = (1. - damping) * cont_freq;

    if (current_time <= xleft) {
        return yleft;
    }
    if (current_time >= xright) {
        return yright;
    }

    else {
        ScalarType normalized_time = (current_time - xleft) / (xright - xleft);

        ScalarType scaled_sigmoid = (sigmoid(sigmoid_param * (normalized_time - 0.5)) - sigmoid(-sigmoid_param * 0.5)) /
                                    (sigmoid(sigmoid_param * 0.5) - sigmoid(-sigmoid_param * 0.5));

        ScalarType smoothed_value = yleft + (yright - yleft) * scaled_sigmoid;

        return smoothed_value;
    }
}

ScalarType sigmoid_smoother_deriv(ScalarType current_time, ScalarType damping_time, ScalarType damping,
                                  ScalarType cont_freq, ScalarType smoother_window, ScalarType sigmoid_param)
{
    ScalarType xleft  = damping_time - smoother_window;
    ScalarType xright = damping_time;

    ScalarType yleft  = cont_freq;
    ScalarType yright = (1. - damping) * cont_freq;

    ScalarType deriv = 0.;

    if (current_time <= xleft || current_time >= xright) {
        deriv = 0.;
    }
    else {
        ScalarType normalized_time       = (current_time - xleft) / (xright - xleft);
        ScalarType normalized_time_deriv = 1. / (xright - xleft);
        deriv = (yright - yleft) / (sigmoid(sigmoid_param * 0.5) - sigmoid(-sigmoid_param * 0.5)) *
                sigmoid(sigmoid_param * (normalized_time - 0.5)) *
                (1. - sigmoid(sigmoid_param * (normalized_time - 0.5))) * sigmoid_param * normalized_time_deriv;
    }

    return deriv;
}

mio::IOResult<void> simulate(ScalarType ide_exponent, size_t finite_difference_order, ScalarType t0_ode,
                             ScalarType t0_ide, ScalarType tmax, ScalarType damping, ScalarType damping_time,
                             ScalarType smoother_window, ScalarType sigmoid_param, std::string save_dir = "",
                             std::string smoother_func_str = "")
{

    using namespace params;
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType dt = pow(10, -ide_exponent);

    // Initial values for SIR model.
    mio::TimeSeries<ScalarType> init_populations((size_t)mio::isir::InfectionState::Count);
    // Initial values for groundtruth.
    mio::TimeSeries<ScalarType> groundtruth_ts((size_t)mio::isir::InfectionState::Count);

    mio::UncertainContactMatrix<ScalarType> contact_matrix = scale_contact_matrix(damping, damping_time);
    Vec vec_init(Vec::Constant((size_t)mio::isir::InfectionState::Count, 0.));

    if (smoother_func_str == "smoothercos") {
        vec_init[(size_t)mio::isir::InfectionState::Susceptible] =
            contact_matrix.get_cont_freq_mat().get_matrix_at(SimulationTime<ScalarType>(ScalarType(0.) * dt))(0, 0);
        vec_init[(size_t)mio::isir::InfectionState::Infected] =
            smoothercos_deriv(0., contact_matrix, damping_time, smoother_window);
        vec_init[(size_t)mio::isir::InfectionState::Recovered] =
            smoothercos_deriv(0., contact_matrix, damping_time, smoother_window);

        init_populations.add_time_point(t0_ode, vec_init);
        groundtruth_ts.add_time_point(t0_ode, vec_init);

        while (init_populations.get_last_time() < t0_ide - 1e-10) {

            vec_init[(size_t)mio::isir::InfectionState::Susceptible] = contact_matrix.get_cont_freq_mat().get_matrix_at(
                SimulationTime<ScalarType>(init_populations.get_last_time() + dt))(0, 0);
            vec_init[(size_t)mio::isir::InfectionState::Infected] =
                smoothercos_deriv(init_populations.get_last_time() + dt, contact_matrix, damping_time, smoother_window);
            vec_init[(size_t)mio::isir::InfectionState::Recovered] =
                smoothercos_deriv(init_populations.get_last_time() + dt, contact_matrix, damping_time, smoother_window);

            init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
            groundtruth_ts.add_time_point(groundtruth_ts.get_last_time() + dt, vec_init);
        }
    }
    else if (smoother_func_str == "smoothstep") {
        vec_init[(size_t)mio::isir::InfectionState::Susceptible] =
            smoothstep(0., damping_time, damping, cont_freq, smoother_window);
        vec_init[(size_t)mio::isir::InfectionState::Infected] =
            smoothstep_deriv(0., damping_time, damping, cont_freq, smoother_window);
        vec_init[(size_t)mio::isir::InfectionState::Recovered] =
            smoothstep_deriv(0., damping_time, damping, cont_freq, smoother_window);

        init_populations.add_time_point(t0_ode, vec_init);
        groundtruth_ts.add_time_point(t0_ode, vec_init);

        while (init_populations.get_last_time() < t0_ide - 1e-10) {

            vec_init[(size_t)mio::isir::InfectionState::Susceptible] =
                smoothstep(init_populations.get_last_time() + dt, damping_time, damping, cont_freq, smoother_window);
            vec_init[(size_t)mio::isir::InfectionState::Infected] = smoothstep_deriv(
                init_populations.get_last_time() + dt, damping_time, damping, cont_freq, smoother_window);
            vec_init[(size_t)mio::isir::InfectionState::Recovered] = smoothstep_deriv(
                init_populations.get_last_time() + dt, damping_time, damping, cont_freq, smoother_window);

            init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
            groundtruth_ts.add_time_point(groundtruth_ts.get_last_time() + dt, vec_init);
        }
    }
    else if (smoother_func_str == "smoothstep_c2") {
        vec_init[(size_t)mio::isir::InfectionState::Susceptible] =
            smoothstep_c2(0., damping_time, damping, cont_freq, smoother_window);
        vec_init[(size_t)mio::isir::InfectionState::Infected] =
            smoothstep_c2_deriv(0., damping_time, damping, cont_freq, smoother_window);
        vec_init[(size_t)mio::isir::InfectionState::Recovered] =
            smoothstep_c2_deriv(0., damping_time, damping, cont_freq, smoother_window);

        init_populations.add_time_point(t0_ode, vec_init);
        groundtruth_ts.add_time_point(t0_ode, vec_init);

        while (init_populations.get_last_time() < t0_ide - 1e-10) {

            vec_init[(size_t)mio::isir::InfectionState::Susceptible] =
                smoothstep_c2(init_populations.get_last_time() + dt, damping_time, damping, cont_freq, smoother_window);
            vec_init[(size_t)mio::isir::InfectionState::Infected] = smoothstep_c2_deriv(
                init_populations.get_last_time() + dt, damping_time, damping, cont_freq, smoother_window);
            vec_init[(size_t)mio::isir::InfectionState::Recovered] = smoothstep_c2_deriv(
                init_populations.get_last_time() + dt, damping_time, damping, cont_freq, smoother_window);

            init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
            groundtruth_ts.add_time_point(groundtruth_ts.get_last_time() + dt, vec_init);
        }
    }
    else if (smoother_func_str == "sigmoid") {
        vec_init[(size_t)mio::isir::InfectionState::Susceptible] =
            sigmoid_smoother(0., damping_time, damping, cont_freq, smoother_window, sigmoid_param);
        vec_init[(size_t)mio::isir::InfectionState::Infected] =
            sigmoid_smoother_deriv(0., damping_time, damping, cont_freq, smoother_window, sigmoid_param);
        vec_init[(size_t)mio::isir::InfectionState::Recovered] =
            sigmoid_smoother_deriv(0., damping_time, damping, cont_freq, smoother_window, sigmoid_param);

        init_populations.add_time_point(t0_ode, vec_init);
        groundtruth_ts.add_time_point(t0_ode, vec_init);

        while (init_populations.get_last_time() < t0_ide - 1e-10) {

            vec_init[(size_t)mio::isir::InfectionState::Susceptible] =
                sigmoid_smoother(init_populations.get_last_time() + dt, damping_time, damping, cont_freq,
                                 smoother_window, sigmoid_param);
            vec_init[(size_t)mio::isir::InfectionState::Infected] =
                sigmoid_smoother_deriv(init_populations.get_last_time() + dt, damping_time, damping, cont_freq,
                                       smoother_window, sigmoid_param);
            vec_init[(size_t)mio::isir::InfectionState::Recovered] =
                sigmoid_smoother_deriv(init_populations.get_last_time() + dt, damping_time, damping, cont_freq,
                                       smoother_window, sigmoid_param);

            init_populations.add_time_point(init_populations.get_last_time() + dt, vec_init);
            groundtruth_ts.add_time_point(groundtruth_ts.get_last_time() + dt, vec_init);
        }
    }

    // Initialize model.
    mio::isir::ModelSmootherCos model(std::move(init_populations), std::move(groundtruth_ts), finite_difference_order,
                                      damping_time, damping, cont_freq, smoother_window, sigmoid_param);

    model.parameters.get<mio::isir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    // Carry out simulation.
    mio::isir::SimulationSmootherCos sim(model, dt);
    sim.advance(tmax, smoother_func_str);

    if (!save_dir.empty()) {
        // Save compartments.
        mio::TimeSeries<ScalarType> compartments = sim.get_result();
        auto save_result_status_ide =
            mio::save_result({compartments}, {0}, num_agegroups,
                             save_dir + "result_ide_dt=1e-" + fmt::format("{:.0f}", ide_exponent) + ".h5");

        mio::TimeSeries<ScalarType> groundtruth = sim.get_groundtruth();
        auto save_result_status_groundtruth =
            mio::save_result({groundtruth}, {0}, num_agegroups,
                             save_dir + "groundtruth_dt=1e-" + fmt::format("{:.0f}", ide_exponent) + ".h5");
    }

    return mio::success();
}

int main()
{
    using namespace params;

    ScalarType t0_ode = 0.;
    ScalarType t0_ide = 5.;
    ScalarType tmax   = 20.;

    ScalarType damping      = 0.2;
    ScalarType damping_time = 15.;

    size_t finite_difference_order = 4;

    ScalarType smoother_window = 5.;

    ScalarType sigmoid_param = 5.;

    // std::string smoother_func_str = "smoothercos";
    // std::string smoother_func_str = "smoothstep";
    // std::string smoother_func_str = "smoothstep_c2";
    std::string smoother_func_str = "sigmoid";

    std::vector<ScalarType> exponents = {0., 1., 2., 3., 4.};

    std::string save_dir =
        fmt::format("../../simulation_results/2026-05-05/{}_fdordercontacts={}_smootherwindow={}_sigmoidparam={}/"
                    "nonconst_contacts_t0={}_tinit={}_tmax={}_dampingtime={}_damping={}_contfreq={}/",
                    smoother_func_str, finite_difference_order, smoother_window, sigmoid_param, t0_ode, t0_ide, tmax,
                    damping_time, damping, cont_freq);

    // Make folder if not existent yet.
    boost::filesystem::path dir(save_dir);
    boost::filesystem::create_directories(dir);

    // Do IDE simulations.

    for (ScalarType exponent : exponents) {
        std::cout << std::endl;
        mio::IOResult<void> result =
            simulate(exponent, finite_difference_order, t0_ode, t0_ide, tmax, damping, damping_time, smoother_window,
                     sigmoid_param, save_dir, smoother_func_str);
    }
}