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
                                   ScalarType smoother_window, ScalarType k)
    : parameters{Parameters()}
    , populations{std::move(populations_init)}
    , groundtruth(std::move(groundtruth_ts))
    , m_finite_difference_order(fd_order)
    , m_damping_time(damping_time)
    , m_damping(damping)
    , m_cont_freq(cont_freq)
    , m_smoother_window(smoother_window)
    , m_k(k)
{
    assert(m_finite_difference_order > 0);
}

ScalarType ModelSmootherCos::compute_deriv_numerical(ScalarType dt, ScalarType current_time,
                                                     std::function<ScalarType(ScalarType)> smoother_func)
{
    ScalarType deriv = 0;

    if (m_finite_difference_order == 1) {
        deriv = (smoother_func(current_time) - smoother_func(current_time - dt)) / dt;
    }

    if (m_finite_difference_order == 2) {
        deriv = (3 * smoother_func(current_time) - 4 * smoother_func(current_time - dt) +
                 1 * smoother_func(current_time - 2. * dt)) /
                (2 * dt);
    }

    if (m_finite_difference_order == 3) {
        deriv = (11 * smoother_func(current_time) - 18 * smoother_func(current_time - dt) +
                 9 * smoother_func(current_time - 2 * dt) - 2 * smoother_func(current_time - 3 * dt)) /
                (6 * dt);
    }

    if (m_finite_difference_order == 4) {
        deriv = (25. * smoother_func(current_time) - 48. * smoother_func(current_time - dt) +
                 36. * smoother_func(current_time - 2. * dt) - 16. * smoother_func(current_time - 3. * dt) +
                 3. * smoother_func(current_time - 4. * dt)) /
                (12. * dt);
    }

    if (m_finite_difference_order == 5) {
        deriv = (137. * smoother_func(current_time) - 300. * smoother_func(current_time - dt) +
                 300. * smoother_func(current_time - 2. * dt) - 200. * smoother_func(current_time - 3. * dt) +
                 75. * smoother_func(current_time - 4. * dt) - 12. * smoother_func(current_time - 5. * dt)) /
                (60. * dt);
    }

    return deriv;
}

ScalarType ModelSmootherCos::smoothercos_via_contacts(ScalarType current_time)
{
    return parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
        SimulationTime<ScalarType>(current_time))(0, 0);
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

void ModelSmootherCos::approximate_smoothercos(ScalarType dt, ScalarType current_time)
{

    // S analytically in Susceptible compartment via contact matrix.
    populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = smoothercos_via_contacts(current_time);

    // S' analytically in Infected compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Infected] = smoothercos_deriv(current_time);

    // S numerically in Recovered compartment via SmootherCOsine in contact_matrix.
    populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] =
        compute_deriv_numerical(dt, current_time, [this](ScalarType t) {
            return smoothercos_via_contacts(t);
        });
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

void ModelSmootherCos::approximate_smoothstep(ScalarType dt, ScalarType current_time)
{
    // S analytically in Susceptible compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = smoothstep(current_time);

    // S' analytically in Infected compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Infected] = smoothstep_deriv(current_time);

    // S numerically in Recovered compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] =
        compute_deriv_numerical(dt, current_time, [this](ScalarType t) {
            return smoothstep(t);
        });
}

ScalarType ModelSmootherCos::smoothstep_c2(ScalarType current_time)
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
            yleft + (yright - yleft) * (6. * std::pow(normalized_time, 5) - 15. * std::pow(normalized_time, 4) +
                                        10. * std::pow(normalized_time, 3));

        return smoothed_value;
    }
}

ScalarType ModelSmootherCos::smoothstep_c2_deriv(ScalarType current_time)
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
        deriv                            = (yright - yleft) *
                (30. * std::pow(normalized_time, 4) - 60. * std::pow(normalized_time, 3) +
                 30 * std::pow(normalized_time, 2)) *
                normalized_time_deriv;
    }

    return deriv;
}

void ModelSmootherCos::approximate_smoothstep_c2(ScalarType dt, ScalarType current_time)
{
    // S analytically in Susceptible compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = smoothstep_c2(current_time);

    // S' analytically in Infected compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Infected] = smoothstep_c2_deriv(current_time);

    // S numerically in Recovered compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] =
        compute_deriv_numerical(dt, current_time, [this](ScalarType t) {
            return smoothstep_c2(t);
        });
}

ScalarType ModelSmootherCos::smoothstep_c4(ScalarType current_time)
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
            yleft + (yright - yleft) * (126. * std::pow(normalized_time, 5) - 420. * std::pow(normalized_time, 6) +
                                        540. * std::pow(normalized_time, 7) - 315. * std::pow(normalized_time, 8) +
                                        70. * std::pow(normalized_time, 9));

        return smoothed_value;
    }
}

ScalarType ModelSmootherCos::smoothstep_c4_deriv(ScalarType current_time)
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
        deriv                            = (yright - yleft) *
                (630. * std::pow(normalized_time, 4) - 2520. * std::pow(normalized_time, 5) +
                 3780. * std::pow(normalized_time, 6) - 2520. * std::pow(normalized_time, 7) +
                 630. * std::pow(normalized_time, 8)) *
                normalized_time_deriv;
    }

    return deriv;
}

void ModelSmootherCos::approximate_smoothstep_c4(ScalarType dt, ScalarType current_time)
{
    // S analytically in Susceptible compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = smoothstep_c4(current_time);

    // S' analytically in Infected compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Infected] = smoothstep_c4_deriv(current_time);

    // S numerically in Recovered compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] =
        compute_deriv_numerical(dt, current_time, [this](ScalarType t) {
            return smoothstep_c4(t);
        });
}

ScalarType ModelSmootherCos::sigmoid(ScalarType x)
{
    return 1. / (1 + std::exp(-x));
}

ScalarType ModelSmootherCos::sigmoid_smoother(ScalarType current_time)
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

        ScalarType scaled_sigmoid =
            (sigmoid(m_k * (normalized_time - 0.5)) - sigmoid(-m_k * 0.5)) / (sigmoid(m_k * 0.5) - sigmoid(-m_k * 0.5));

        ScalarType smoothed_value = yleft + (yright - yleft) * scaled_sigmoid;

        return smoothed_value;
    }
}

ScalarType ModelSmootherCos::sigmoid_smoother_deriv(ScalarType current_time)
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
        deriv = (yright - yleft) / (sigmoid(m_k * 0.5) - sigmoid(-m_k * 0.5)) * sigmoid(m_k * (normalized_time - 0.5)) *
                (1. - sigmoid(m_k * (normalized_time - 0.5))) * m_k * normalized_time_deriv;
    }

    return deriv;
}

void ModelSmootherCos::approximate_sigmoid_smoother(ScalarType dt, ScalarType current_time)
{
    // S analytically in Susceptible compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = sigmoid_smoother(current_time);

    // S' analytically in Infected compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Infected] = sigmoid_smoother_deriv(current_time);

    // S numerically in Recovered compartment.
    populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] =
        compute_deriv_numerical(dt, current_time, [this](ScalarType t) {
            return sigmoid_smoother(t);
        });
}

void ModelSmootherCos::set_groundtruth(ScalarType current_time, std::string smoother_func_str)
{
    if (smoother_func_str == "smoothercos") {
        // Set to smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = smoothercos(current_time);
        // Set to derivative of smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Infected]  = smoothercos_deriv(current_time);
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Recovered] = smoothercos_deriv(current_time);
    }

    else if (smoother_func_str == "smoothstep") {
        // Set to smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = smoothstep(current_time);
        // Set to derivative of smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Infected]  = smoothstep_deriv(current_time);
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Recovered] = smoothstep_deriv(current_time);
    }
    else if (smoother_func_str == "smoothstep_c2") {
        // Set to smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = smoothstep_c2(current_time);
        // Set to derivative of smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Infected]  = smoothstep_c2_deriv(current_time);
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Recovered] = smoothstep_c2_deriv(current_time);
    }
    else if (smoother_func_str == "smoothstep_c4") {
        // Set to smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = smoothstep_c4(current_time);
        // Set to derivative of smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Infected]  = smoothstep_c4_deriv(current_time);
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Recovered] = smoothstep_c4_deriv(current_time);
    }
    else if (smoother_func_str == "sigmoid") {
        // Set to smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = sigmoid_smoother(current_time);
        // Set to derivative of smoother function.
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Infected]  = sigmoid_smoother_deriv(current_time);
        groundtruth.get_last_value()[(Eigen::Index)InfectionState::Recovered] = sigmoid_smoother_deriv(current_time);
    }
}

} // namespace isir
} // namespace mio
