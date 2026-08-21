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

#include "ide_sir/infection_state.h"
#include "ide_sir/model.h"
#include "ide_sir/simulation.h"
#include "ide_sir/parameters.h"
#include "memilio/config.h"
#include "memilio/epidemiology/damping.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/time_series.h"
#include <Eigen/src/Core/util/Meta.h>
#include <cmath>

namespace mio
{
namespace isir
{

using Vec = mio::TimeSeries<ScalarType>::Vector;

void SimulationMessinaExtendedDetailedInit::advance(ScalarType tmax, bool kahan, bool more_precise_s_deriv,
                                                    ScalarType cutoff_window, size_t fd_order_contacts,
                                                    ScalarType damping_time, ScalarType smoother_window,
                                                    size_t smoothstep_order)
{
    if ((smoothstep_order != 0) && (smoothstep_order != 1) && (smoothstep_order != 3) && (smoothstep_order != 4)) {
        throw std::invalid_argument("smoothstep_order must be 0, 3, or 4");
    }
    // Get index of t0, i.e. index of last time point of given initial values.
    size_t t0_index = m_model->populations.get_num_time_points() - 1;

    // Set vector with values of transition distribution and parameters, respectively.
    m_model->set_transitiondistribution_vector(m_dt, tmax, t0_index);
    m_model->set_parameter_vectors(m_dt, tmax, t0_index);

    mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {}.", m_model->populations.get_last_time(),
                  tmax, m_dt);

    // Compute S' for t_0,..., t_{n0-1}.
    // We set S'(0) due to lack of knowledge of previous values of S.
    // The corresponding flow is then given by -S'.
    // TODO: Initialize S'(0) in a different way?

    if (m_model->flows.get_num_time_points() == 0) {

        // ScalarType first_flow_approx = m_model->parameters.get<TransmissionProbabilityOnContact>().eval(0) *
        //                                m_model->parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
        //                                    SimulationTime<ScalarType>(0.))(0, 0) *
        //                                (m_model->populations.get_value(0)[(Eigen::Index)InfectionState::Infected] +
        //                                 m_model->populations.get_value(0)[(Eigen::Index)InfectionState::Recovered]) /
        //                                m_model->populations.get_value(0).sum();

        ScalarType first_flow_approx = m_model->populations.get_value(0)[(Eigen::Index)InfectionState::Infected];

        ScalarType t_flows_init = m_model->populations.get_time(0);
        m_model->flows.add_time_point(t_flows_init, TimeSeries<ScalarType>::Vector::Constant(
                                                        (size_t)InfectionTransition::Count, first_flow_approx));
        std::cout << "Flows first tp: " << m_model->flows.get_time(0) << std::endl;
        // Compute S'(t) for t_1,..., t_{n0-1} with backwards difference operator. The corresponding flow is then given by -S'.
        for (size_t i = 1; i < (size_t)m_model->populations.get_num_time_points(); i++) {

            ScalarType increment         = m_dt - m_summation_error_flows_init;
            ScalarType t_temp            = t_flows_init + increment;
            m_summation_error_flows_init = (t_temp - t_flows_init) - increment;
            // std::cout << t << ", " << m_summation_error << std::endl;
            t_flows_init = t_temp;

            m_model->flows.add_time_point(
                t_flows_init, TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));

            m_model->compute_S_deriv(m_div_dt, i);
        }

        if (more_precise_s_deriv) {
            // Only the extra buffer at the front was needed to get a more accurate finite difference approximation
            // of S' close to t_init. Drop it now from both flows and populations so that the two time series stay
            // aligned index-by-index (compute_S_deriv() and compute_I_and_R() rely on this alignment when deriving
            // their time_point_index from the number of time points already present in flows).
            size_t num_buffer_points = std::round(cutoff_window / m_dt);
            for (size_t i = 0; i < num_buffer_points; i++) {
                m_model->flows.remove_time_point(0);
                m_model->populations.remove_time_point(0);
            }
            t0_index -= num_buffer_points;
        }
        std::cout << "flows first tp: " << m_model->flows.get_time(0) << std::endl;
    }

    ScalarType t_pop = m_model->populations.get_last_time();

    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        ScalarType increment  = m_dt - m_summation_error_pop;
        ScalarType t_temp     = t_pop + increment;
        m_summation_error_pop = (t_temp - t_pop) - increment;
        // std::cout << t << ", " << m_summation_error << std::endl;
        t_pop = t_temp;

        // Add new time point to populations.
        m_model->populations.add_time_point(t_pop, Vec::Constant((size_t)InfectionState::Count, 0.));

        // Print time.
        if (floating_point_equal(std::remainder(10 * m_model->populations.get_last_time(), tmax), 0., 1e-7)) {
            std::cout << "Time pop: " << m_model->populations.get_last_time() << std::endl;
        }

        // std::cout << "Time pop: " << m_model->populations.get_last_time() << std::endl;

        // Compute Susceptibles.
        size_t num_time_points = m_model->populations.get_num_time_points();

        size_t num_iterations =
            m_model->compute_S(m_model->populations.get_value(num_time_points - 2)[(size_t)InfectionState::Susceptible],
                               m_dt, fd_order_contacts, kahan, damping_time, smoother_window, smoothstep_order);

        if (num_iterations > m_max_number_iterations) {
            m_max_number_iterations = num_iterations;
        }
    }

    // Compute S' as well as I and R.
    ScalarType t_flows = m_model->flows.get_last_time();
    while (m_model->flows.get_last_time() < tmax - 1e-10) {

        ScalarType increment    = m_dt - m_summation_error_flows;
        ScalarType t_temp       = t_flows + increment;
        m_summation_error_flows = (t_temp - t_flows) - increment;
        // std::cout << t << ", " << m_summation_error << std::endl;
        t_flows = t_temp;

        m_model->flows.add_time_point(t_flows, Vec::Constant((size_t)InfectionTransition::Count, 0.));

        if (floating_point_equal(std::remainder(10 * m_model->flows.get_last_time(), tmax), 0., 1e-7)) {
            std::cout << "Time flows: " << m_model->flows.get_last_time() << std::endl;
        }

        // Compute S'.
        m_model->compute_S_deriv(m_div_dt);

        // Compute I and R.
        m_model->compute_I_and_R(m_dt, kahan);
    }

    std::cout << "SIR: " << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] << std::endl;

    // std::cout << "Total population at start of simulation is "
    //           << m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Susceptible] +
    //                  m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Infected] +
    //                  m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Recovered]
    //           << " and at end "
    //           << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] +
    //                  m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] +
    //                  m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered]
    //           << std::endl;

    std::cout << "Difference in total population: "
              << m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Susceptible] +
                     m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Infected] +
                     m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Recovered] -
                     (m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] +
                      m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] +
                      m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered])
              << std::endl;

    std::cout << "Max number of iterations throughout simulation was " << m_max_number_iterations << std::endl;

    // std::cout << "t0_index: " << t0_index << std::endl;

    auto file = m_model->populations.export_csv("populations_ide.csv");
    std::cout << std::endl;
}

void SimulationMessinaExtendedDetailedInit::advance_S_deriv_analytical(ScalarType tmax, size_t fd_order_contacts,
                                                                       ScalarType damping_time)
{
    // Get index of t0, i.e. index of last time point of given initial values.
    size_t t0_index = m_model->populations.get_num_time_points() - 1;

    // Set vector with values of transition distribution and parameters, respectively.
    m_model->set_transitiondistribution_vector(m_dt, tmax, t0_index);
    m_model->set_parameter_vectors(m_dt, tmax, t0_index);

    mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {}.", m_model->populations.get_last_time(),
                  tmax, m_dt);

    // Compute S' for t_0,..., t_{n0-1}.
    // We set S'(0) due to lack of knowledge of previous values of S.
    // The corresponding flow is then given by -S'.
    // TODO: Initialize S'(0) in a different way?

    if (m_model->flows.get_num_time_points() == 0) {

        ScalarType first_flow_approx = m_model->populations.get_value(0)[(Eigen::Index)InfectionState::Infected];

        ScalarType t_flows_init = m_model->populations.get_time(0);
        m_model->flows.add_time_point(t_flows_init, TimeSeries<ScalarType>::Vector::Constant(
                                                        (size_t)InfectionTransition::Count, first_flow_approx));
        std::cout << "Flows first tp: " << m_model->flows.get_time(0) << std::endl;
        // Compute S'(t) for t_1,..., t_{n0-1} with backwards difference operator. The corresponding flow is then given by -S'.
        for (size_t i = 1; i < (size_t)m_model->populations.get_num_time_points(); i++) {

            ScalarType increment         = m_dt - m_summation_error_flows_init;
            ScalarType t_temp            = t_flows_init + increment;
            m_summation_error_flows_init = (t_temp - t_flows_init) - increment;
            // std::cout << t << ", " << m_summation_error << std::endl;
            t_flows_init = t_temp;

            m_model->flows.add_time_point(
                t_flows_init, TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));

            m_model->compute_S_deriv_analytical(m_dt, i);
        }
    }

    ScalarType t_pop = m_model->populations.get_last_time();

    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        ScalarType increment  = m_dt - m_summation_error_pop;
        ScalarType t_temp     = t_pop + increment;
        m_summation_error_pop = (t_temp - t_pop) - increment;
        // std::cout << t << ", " << m_summation_error << std::endl;
        t_pop = t_temp;

        // Add new time point to populations.
        m_model->populations.add_time_point(t_pop, Vec::Constant((size_t)InfectionState::Count, 0.));

        // Print time.
        if (floating_point_equal(std::remainder(10 * m_model->populations.get_last_time(), tmax), 0., 1e-7)) {
            std::cout << "Time pop: " << m_model->populations.get_last_time() << std::endl;
        }

        // Compute Susceptibles.
        size_t num_time_points = m_model->populations.get_num_time_points();

        size_t num_iterations =
            m_model->compute_S(m_model->populations.get_value(num_time_points - 2)[(size_t)InfectionState::Susceptible],
                               m_dt, fd_order_contacts, damping_time);

        if (num_iterations > m_max_number_iterations) {
            m_max_number_iterations = num_iterations;
        }
    }

    // Compute S' as well as I and R.
    ScalarType t_flows = m_model->flows.get_last_time();
    while (m_model->flows.get_last_time() < tmax - 1e-10) {

        long double increment   = m_dt - m_summation_error_flows;
        long double t_temp      = t_flows + increment;
        m_summation_error_flows = (t_temp - (long double)t_flows) - increment;
        // std::cout << t << ", " << m_summation_error << std::endl;
        t_flows = (ScalarType)t_temp;

        m_model->flows.add_time_point(t_flows, Vec::Constant((size_t)InfectionTransition::Count, 0.));

        if (floating_point_equal(std::remainder(10 * m_model->flows.get_last_time(), tmax), 0., 1e-7)) {
            std::cout << "Time flows: " << m_model->flows.get_last_time() << std::endl;
        }

        // Compute S'.
        int last_tp_index = m_model->flows.get_num_time_points() - 1;
        m_model->compute_S_deriv_analytical(m_dt, last_tp_index);

        // Compute I and R.
        m_model->compute_I_and_R(m_dt);
    }

    std::cout << "SIR: " << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] << std::endl;

    std::cout << "Difference in total population: "
              << m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Susceptible] +
                     m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Infected] +
                     m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Recovered] -
                     (m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] +
                      m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] +
                      m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered])
              << std::endl;

    std::cout << "Max number of iterations throughout simulation was " << m_max_number_iterations << std::endl;

    // std::cout << "t0_index: " << t0_index << std::endl;

    auto file = m_model->populations.export_csv("populations_ide.csv");
    std::cout << std::endl;
}

void SimulationMessinaExtendedDetailedInit::advance_S_deriv_fixedpoint(ScalarType tmax, size_t fd_order_contacts,
                                                                       ScalarType damping_time)
{
    // Get index of t0, i.e. index of last time point of given initial values.
    size_t t0_index = m_model->populations.get_num_time_points() - 1;

    // Set vector with values of transition distribution and parameters, respectively.
    m_model->set_transitiondistribution_vector(m_dt, tmax, t0_index);
    m_model->set_parameter_vectors(m_dt, tmax, t0_index);

    mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {}.", m_model->populations.get_last_time(),
                  tmax, m_dt);

    // Compute S' for t_0,..., t_{n0-1}.
    // We set S'(0) due to lack of knowledge of previous values of S.
    // The corresponding flow is then given by -S'.
    // TODO: Initialize S'(0) in a different way?

    m_model->flows.add_time_point(m_model->populations.get_time(0),
                                  TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));
    std::cout << "Flows first tp: " << m_model->flows.get_time(0) << std::endl;
    // Compute S'(t) for t_1,..., t_{n0-1} with backwards difference operator. The corresponding flow is then given by -S'.
    for (size_t i = 1; i < (size_t)m_model->populations.get_num_time_points(); i++) {
        m_model->flows.add_time_point(m_model->flows.get_last_time() + m_dt,
                                      TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));

        m_model->compute_S_deriv(m_dt, i);
    }

    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        // Add new time point to populations.
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));

        // Print time.
        if (floating_point_equal(std::remainder(10 * m_model->populations.get_last_time(), tmax), 0., 1e-7)) {
            std::cout << "Time pop: " << m_model->populations.get_last_time() << std::endl;
        }

        // std::cout << "Time pop: " << m_model->populations.get_last_time() << std::endl;

        // Compute Susceptibles.
        size_t num_time_points = m_model->populations.get_num_time_points();

        size_t num_iterations =
            m_model->compute_S(m_model->populations.get_value(num_time_points - 2)[(size_t)InfectionState::Susceptible],
                               m_dt, fd_order_contacts, damping_time);

        if (num_iterations > m_max_number_iterations) {
            m_max_number_iterations = num_iterations;
        }
    }

    // Compute S' as well as I and R.
    while (m_model->flows.get_last_time() < tmax - 1e-10) {

        m_model->flows.add_time_point(m_model->flows.get_last_time() + m_dt,
                                      Vec::Constant((size_t)InfectionTransition::Count, 0.));

        if (floating_point_equal(std::remainder(10 * m_model->flows.get_last_time(), tmax), 0., 1e-7)) {
            std::cout << "Time flows: " << m_model->flows.get_last_time() << std::endl;
        }

        // Compute S'.
        size_t num_time_points = m_model->flows.get_num_time_points();
        ScalarType s_deriv_init =
            -m_model->flows.get_value(num_time_points - 2)[(size_t)InfectionTransition::SusceptibleToInfected];
        m_model->compute_S_deriv_fixedpoint(s_deriv_init, m_dt);

        // Compute I and R.
        m_model->compute_I_and_R(m_dt);
    }

    std::cout << "SIR: " << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] << std::endl;

    // std::cout << "Total population at start of simulation is "
    //           << m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Susceptible] +
    //                  m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Infected] +
    //                  m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Recovered]
    //           << " and at end "
    //           << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] +
    //                  m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] +
    //                  m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered]
    //           << std::endl;

    std::cout << "Difference in total population: "
              << m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Susceptible] +
                     m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Infected] +
                     m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Recovered] -
                     (m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] +
                      m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] +
                      m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered])
              << std::endl;

    std::cout << "Max number of iterations throughout simulation was " << m_max_number_iterations << std::endl;

    // std::cout << "t0_index: " << t0_index << std::endl;

    auto file = m_model->populations.export_csv("populations_ide.csv");
    std::cout << std::endl;
}

void SimulationMessinaExtendedDetailedInit::advance_reformulated(ScalarType tmax)
{
    // Get index of t0, i.e. index of last time point of given initial values.
    size_t t0_index = m_model->populations.get_num_time_points() - 1;

    // Set vector with values of transition distribution and parameters, respectively.
    m_model->set_transitiondistribution_vector(m_dt, tmax, t0_index);
    m_model->set_parameter_vectors(m_dt, tmax, t0_index);

    mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {}.", m_model->populations.get_last_time(),
                  tmax, m_dt);

    // Compute S' for t_0,..., t_{n0-1}.
    // We set S'(0) due to lack of knowledge of previous values of S.
    // The corresponding flow is then given by -S'.
    // TODO: Initialize S'(0) in a different way?
    m_model->flows.add_time_point(0., TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));

    // Compute S'(t) for t_1,..., t_{n0-1} with backwards difference operator. The corresponding flow is then given by -S'.
    for (size_t i = 1; i < (size_t)m_model->populations.get_num_time_points(); i++) {
        m_model->flows.add_time_point(i * m_dt,
                                      TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));
        m_model->compute_S_deriv(m_dt, i);
    }

    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        // // Print time.
        // if (floating_point_equal(std::remainder(10 * m_model->populations.get_last_time(), tmax), 0., 1e-7)) {
        //     std::cout << "Time pop: " << m_model->populations.get_last_time() << std::endl;
        // }
        std::cout << "Time pop: " << m_model->populations.get_last_time() << std::endl;

        // Add new time point to populations.
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));

        // Compute Susceptibles.
        size_t num_time_points = m_model->populations.get_num_time_points();

        size_t num_iterations = m_model->compute_S_reformulated(
            m_model->populations.get_value(num_time_points - 2)[(size_t)InfectionState::Susceptible], m_dt, t0_index);

        if (num_iterations > m_max_number_iterations) {
            m_max_number_iterations = num_iterations;
        }

        m_model->flows.add_time_point(m_model->flows.get_last_time() + m_dt,
                                      Vec::Constant((size_t)InfectionTransition::Count, 0.));

        // Compute S'.
        m_model->compute_S_deriv(m_dt);

        // Compute I and R.
        m_model->compute_I_and_R(m_dt, m_model->flows.get_num_time_points());
    }

    std::cout << "SIR: " << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] << std::endl;

    std::cout << "Difference in total population: "
              << m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Susceptible] +
                     m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Infected] +
                     m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Recovered] -
                     (m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] +
                      m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] +
                      m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered])
              << std::endl;

    std::cout << "Max number of iterations throughout simulation was " << m_max_number_iterations << std::endl;
}

void SimulationMessinaExtendedDetailedInit::advance_reformulated2(ScalarType tmax, size_t fd_order)
{
    // Get index of t0, i.e. index of last time point of given initial values.
    size_t t0_index = m_model->populations.get_num_time_points() - 1;

    // Set vector with values of transition distribution and parameters, respectively.
    m_model->set_transitiondistribution_vector(m_dt, tmax, t0_index);
    m_model->set_parameter_vectors(m_dt, tmax, t0_index);

    mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {}.", m_model->populations.get_last_time(),
                  tmax, m_dt);

    // Compute S' for t_0,..., t_{n0-1}.
    // We set S'(0) due to lack of knowledge of previous values of S.
    // The corresponding flow is then given by -S'.
    // TODO: Initialize S'(0) in a different way?
    m_model->flows.add_time_point(0., TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));

    // Compute S'(t) for t_1,..., t_{n0-1} with backwards difference operator. The corresponding flow is then given by -S'.
    for (size_t i = 1; i < (size_t)m_model->populations.get_num_time_points(); i++) {
        m_model->flows.add_time_point(i * m_dt,
                                      TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));
        m_model->compute_S_deriv(m_dt, i);
    }

    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        // // Print time.
        // if (floating_point_equal(std::remainder(10 * m_model->populations.get_last_time(), tmax), 0., 1e-7)) {
        //     std::cout << "Time pop: " << m_model->populations.get_last_time() << std::endl;
        // }
        std::cout << "Time pop: " << m_model->populations.get_last_time() << std::endl;

        // Add new time point to populations.
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));

        // Compute Susceptibles.
        size_t num_time_points = m_model->populations.get_num_time_points();

        size_t num_iterations = m_model->compute_S_reformulated2(
            m_model->populations.get_value(num_time_points - 2)[(size_t)InfectionState::Susceptible], m_dt, t0_index,
            fd_order);

        if (num_iterations > m_max_number_iterations) {
            m_max_number_iterations = num_iterations;
        }
    }

    // Compute S' as well as I and R.
    while (m_model->flows.get_last_time() < tmax - 1e-10) {

        if (floating_point_equal(std::remainder(10 * m_model->flows.get_last_time(), tmax), 0., 1e-7)) {
            std::cout << "Time flows: " << m_model->flows.get_last_time() << std::endl;
        }

        m_model->flows.add_time_point(m_model->flows.get_last_time() + m_dt,
                                      Vec::Constant((size_t)InfectionTransition::Count, 0.));

        // Compute S'.
        m_model->compute_S_deriv(m_dt);

        // Compute I and R.
        m_model->compute_I_and_R(m_dt);
    }

    std::cout << "SIR: " << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] << std::endl;

    // std::cout << "Total population at start of simulation is "
    //           << m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Susceptible] +
    //                  m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Infected] +
    //                  m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Recovered]
    //           << " and at end "
    //           << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] +
    //                  m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] +
    //                  m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered]
    //           << std::endl;

    std::cout << "Difference in total population: "
              << m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Susceptible] +
                     m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Infected] +
                     m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Recovered] -
                     (m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] +
                      m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] +
                      m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered])
              << std::endl;

    std::cout << "Max number of iterations throughout simulation was " << m_max_number_iterations << std::endl;

    // std::cout << "t0_index: " << t0_index << std::endl;

    auto file = m_model->populations.export_csv("populations_ide.csv");
    std::cout << std::endl;
}

} // namespace isir
} // namespace mio
