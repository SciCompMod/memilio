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
#include "ide_sir/parameters.h"
#include "ide_sir/infection_state.h"
#include "ide_sir/gregory_weights.h"
#include "memilio/config.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/time_series.h"
#include <Eigen/src/Core/util/Meta.h>
#include <boost/numeric/odeint/util/ublas_wrapper.hpp>
#include <cmath>

namespace mio
{
namespace isir
{
ModelSmootherCos::ModelSmootherCos(TimeSeries<ScalarType>&& populations_init, TimeSeries<ScalarType>&& groundtruth_ts,
                                   size_t fd_order, ScalarType damping_time, ScalarType damping, ScalarType cont_freq,
                                   ScalarType smoother_window)
    : parameters{Parameters()}
    , populations{std::move(populations_init)}
    , groundtruth(std::move(groundtruth_ts))
    , m_finite_difference_order(fd_order)
    , m_damping_time(damping_time)
    , m_damping(damping)
    , m_cont_freq(cont_freq)
    , m_smoother_window(smoother_window)
{
    assert(m_gregory_order > 0);
    assert(m_finite_difference_order > 0);
}

ScalarType ModelSmootherCos::smoothercos(ScalarType current_time)
{
    ScalarType xleft  = m_damping_time - m_smoother_window;
    ScalarType xright = m_damping_time;

    ScalarType yleft  = m_cont_freq;
    ScalarType yright = (1. - m_damping) * m_cont_freq;

    if (current_time <= xleft) {
        return yleft;
    }

    if (current_time >= xright) {
        return yright;
    }

    else {
        return 0.5 * (yleft - yright) *
                   std::cos(std::numbers::pi_v<ScalarType> / (xright - xleft) * (current_time - xleft)) +
               0.5 * (yleft + yright);
    }
}

ScalarType ModelSmootherCos::smoothercos_deriv(ScalarType current_time)
{
    ScalarType xleft  = m_damping_time - m_smoother_window;
    ScalarType xright = m_damping_time;

    if (current_time <= xleft || current_time >= xright) {
        return 0.;
    }

    ScalarType yleft =
        parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(SimulationTime<ScalarType>(xleft))(0, 0);
    ScalarType yright =
        parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(SimulationTime<ScalarType>(xright))(0, 0);

    ScalarType deriv = -0.5 * (yleft - yright) / (xright - xleft) * std::numbers::pi_v<ScalarType> *
                       std::sin(std::numbers::pi_v<ScalarType> / (xright - xleft) * (current_time - xleft));

    return deriv;
}

// void ModelSmootherCos::compute_S()
// {
//     ScalarType current_time = populations.get_last_time();
//     populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] =
//         parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(SimulationTime<ScalarType>(current_time))(
//             0, 0);
// }

// void ModelSmootherCos::compute_S_deriv_analytical(ScalarType current_time)
// {

//     populations.get_last_value()[(Eigen::Index)InfectionState::Infected] = smoothercos_deriv(current_time);
// }

void ModelSmootherCos::smoothercos_deriv_numerical(ScalarType dt, size_t j)
{
    ScalarType deriv = 0;

    if (m_finite_difference_order == 1) {
        deriv = (parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                     SimulationTime<ScalarType>(ScalarType(j) * dt))(0, 0) -
                 parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                     SimulationTime<ScalarType>((ScalarType(j) - 1.) * dt))(0, 0)) /
                dt;
    }

    if (m_finite_difference_order == 2) {
        deriv = (3 * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                         SimulationTime<ScalarType>(ScalarType(j) * dt))(0, 0) -
                 4 * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                         SimulationTime<ScalarType>((ScalarType(j) - 1.) * dt))(0, 0) +
                 1 * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                         SimulationTime<ScalarType>(((ScalarType)j - 2.) * dt))(0, 0)) /
                (2 * dt);
    }

    if (m_finite_difference_order == 3) {
        deriv = (11 * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                          SimulationTime<ScalarType>(ScalarType(j) * dt))(0, 0) -
                 18 * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                          SimulationTime<ScalarType>((ScalarType(j) - 1.) * dt))(0, 0) +
                 9 * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                         SimulationTime<ScalarType>(((ScalarType)j - 2.) * dt))(0, 0) -
                 2 * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                         SimulationTime<ScalarType>(((ScalarType)j - 3.) * dt))(0, 0)) /
                (6 * dt);
    }

    if (m_finite_difference_order == 4) {
        deriv = (25. * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                           SimulationTime<ScalarType>((ScalarType)j * dt))(0, 0) -
                 48. * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                           SimulationTime<ScalarType>(((ScalarType)j - 1.) * dt))(0, 0) +
                 36. * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                           SimulationTime<ScalarType>(((ScalarType)j - 2.) * dt))(0, 0) -
                 16. * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                           SimulationTime<ScalarType>(((ScalarType)j - 3.) * dt))(0, 0) +
                 3. * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                          SimulationTime<ScalarType>(((ScalarType)j - 4.) * dt))(0, 0)) /
                (12. * dt);
    }

    if (m_finite_difference_order == 5) {
        deriv = (137. * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                            SimulationTime<ScalarType>((ScalarType)j * dt))(0, 0) -
                 300. * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                            SimulationTime<ScalarType>(((ScalarType)j - 1.) * dt))(0, 0) +
                 300. * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                            SimulationTime<ScalarType>(((ScalarType)j - 2.) * dt))(0, 0) -
                 200. * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                            SimulationTime<ScalarType>(((ScalarType)j - 3.) * dt))(0, 0) +
                 75. * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                           SimulationTime<ScalarType>(((ScalarType)j - 4.) * dt))(0, 0) -
                 12. * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                           SimulationTime<ScalarType>(((ScalarType)j - 5.) * dt))(0, 0)) /
                (60. * dt);
    }

    populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] = deriv;
}

void ModelSmootherCos::smoothercos_deriv_numerical(ScalarType dt)
{
    // Use the number of time points to determine time_point_index, hence we are calculating S deriv for last time point.
    size_t j = populations.get_num_time_points() - 1;
    smoothercos_deriv_numerical(dt, j);
}

void ModelSmootherCos::approximate_smoothercos(ScalarType dt, ScalarType current_time)
{

    // S analytically in Susceptible compartment via contact matrix.
    populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] =
        parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(SimulationTime<ScalarType>(current_time))(
            0, 0);

    // S' analytically in Infected compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Infected] = smoothercos_deriv(current_time);

    // S numerically in Recovered compartment.
    smoothercos_deriv_numerical(dt);
}

ScalarType ModelSmootherCos::smoothstep(ScalarType current_time)
{

    ScalarType xleft  = m_damping_time - m_smoother_window;
    ScalarType xright = m_damping_time;

    ScalarType yleft  = m_cont_freq;
    ScalarType yright = (1. - m_damping) * m_cont_freq;

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

ScalarType ModelSmootherCos::smoothstep_deriv(ScalarType current_time)
{
    ScalarType xleft  = m_damping_time - m_smoother_window;
    ScalarType xright = m_damping_time;

    ScalarType yleft  = m_cont_freq;
    ScalarType yright = (1. - m_damping) * m_cont_freq;

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

void ModelSmootherCos::smoothstep_deriv_numerical(ScalarType current_time, ScalarType dt)
{
    ScalarType deriv = 0;

    if (m_finite_difference_order == 1) {
        deriv = (smoothstep(current_time) - smoothstep(current_time - dt)) / dt;
    }

    if (m_finite_difference_order == 2) {
        deriv = (3 * smoothstep(current_time) - 4 * smoothstep(current_time - dt) +
                 1 * smoothstep(current_time - 2. * dt)) /
                (2 * dt);
    }

    if (m_finite_difference_order == 3) {
        deriv = (11 * smoothstep(current_time) - 18 * smoothstep(current_time - dt) +
                 9 * smoothstep(current_time - 2 * dt) - 2 * smoothstep(current_time - 3 * dt)) /
                (6 * dt);
    }

    if (m_finite_difference_order == 4) {
        deriv = (25. * smoothstep(current_time) - 48. * smoothstep(current_time - dt) +
                 36. * smoothstep(current_time - 2. * dt) - 16. * smoothstep(current_time - 3. * dt) +
                 3. * smoothstep(current_time - 4. * dt)) /
                (12. * dt);
    }

    if (m_finite_difference_order == 5) {
        deriv = (137. * smoothstep(current_time) - 300. * smoothstep(current_time - dt) +
                 300. * smoothstep(current_time - 2. * dt) - 200. * smoothstep(current_time - 3. * dt) +
                 75. * smoothstep(current_time - 4. * dt) - 12. * smoothstep(current_time - 5. * dt)) /
                (60. * dt);
    }

    populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] = deriv;
}

void ModelSmootherCos::approximate_smoothstep(ScalarType dt, ScalarType current_time)
{
    // S analytically in Susceptible compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = smoothstep(current_time);

    // S' analytically in Infected compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Infected] = smoothstep_deriv(current_time);

    // S numerically in Recovered compartment.
    smoothstep_deriv_numerical(current_time, dt);
}

void ModelSmootherCos::set_groundtruth(ScalarType current_time, bool smoothercos_func)
{
    if (smoothercos_func) {
        // Set to smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = smoothercos(current_time);
        // Set to derivative of smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Infected]  = smoothercos_deriv(current_time);
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Recovered] = smoothercos_deriv(current_time);
    }

    else {
        // Set to smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = smoothstep(current_time);
        // Set to derivative of smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Infected]  = smoothstep_deriv(current_time);
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Recovered] = smoothstep_deriv(current_time);
    }
}

} // namespace isir
} // namespace mio
