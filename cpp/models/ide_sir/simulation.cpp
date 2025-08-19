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
#include "memilio/config.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/time_series.h"
#include <cmath>

namespace mio
{
namespace isir
{

using Vec = mio::TimeSeries<ScalarType>::Vector;

// In this advance function, we only consider S as in the paper by Messina et al.
void SimulationMessina::advance_messina(ScalarType tmax)
{
    m_model->set_transitiondistribution_vector(m_dt, tmax);
    m_model->set_parameter_vectors(m_dt, tmax);

    mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {}.", m_model->populations.get_last_time(),
                  tmax, m_dt);

    while (m_model->populations.get_last_time() < tmax - 1e-10) {
        if (floating_point_equal(std::remainder(10 * m_model->populations.get_last_time(), tmax), 0., 1e-7)) {
            std::cout << "Time: " << m_model->populations.get_last_time() << std::endl;
        }

        // Add new time point to populations.
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));

        // Compute Susceptibles.
        size_t num_time_points = m_model->populations.get_num_time_points();
        // std::cout << "S_init "
        //           << m_model->populations.get_value(num_time_points - 2)[(size_t)InfectionState::Susceptible]
        //           << " at time " << m_model->populations.get_time(num_time_points - 2) << std::endl;
        size_t num_iterations = m_model->compute_S(
            m_model->populations.get_value(num_time_points - 2)[(size_t)InfectionState::Susceptible], m_dt);
        // m_model->populations.print_table();

        if (num_iterations > m_max_number_iterations) {
            m_max_number_iterations = num_iterations;
        }
    }
    std::cout << "Max number of iterations throughout simulation was " << m_max_number_iterations << std::endl;
}

// In this advance function, we only consider S as in the paper by Messina et al.
void SimulationMessinaExtended::advance_messina(ScalarType tmax)
{

    size_t t0_index = m_model->populations.get_num_time_points() - 1;
    std::cout << "t0_index: " << t0_index << std::endl;

    m_model->set_transitiondistribution_vector(m_dt, tmax);
    m_model->set_parameter_vectors(m_dt, tmax);

    mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {}.", m_model->populations.get_last_time(),
                  tmax, m_dt);

    // Compute S' for t_0,..., t_{n0-1}.
    // Compute S'(0) with forward difference operator due to lack of knowledge of previous values of S.
    // The corresponding flow is then given by -S'.
    m_model->flows.add_time_point(0., TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));
    m_model->flows.get_value(0)[(Eigen::Index)InfectionTransition::SusceptibleToInfected] = 0.;

    // Compute S'(t) for t_1,..., t_{n0-1} with backwards difference operator. The corresponding flow is then given by -S'.
    for (size_t i = 1; i < (size_t)m_model->populations.get_num_time_points(); i++) {
        m_model->flows.add_time_point(i * m_dt,
                                      TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));
        m_model->flows.get_value(i)[(Eigen::Index)InfectionTransition::SusceptibleToInfected] =
            -(m_model->populations.get_value(i)[(Eigen::Index)InfectionState::Susceptible] -
              m_model->populations.get_value(i - 1)[(Eigen::Index)InfectionState::Susceptible]) /
            m_dt;
    }

    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        // Add new time point to populations.
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));
        m_model->flows.add_time_point(m_model->flows.get_last_time() + m_dt,
                                      Vec::Constant((size_t)InfectionTransition::Count, 0.));

        // Print time.
        if (floating_point_equal(std::remainder(10 * m_model->populations.get_last_time(), tmax), 0., 1e-7)) {
            std::cout << "Time: " << m_model->populations.get_last_time() << std::endl;
        }

        // Compute Susceptibles.
        size_t num_time_points = m_model->populations.get_num_time_points();

        size_t num_iterations = m_model->compute_S(
            m_model->populations.get_value(num_time_points - 2)[(size_t)InfectionState::Susceptible], m_dt);

        if (num_iterations > m_max_number_iterations) {
            m_max_number_iterations = num_iterations;
        }

        // Compute S'.
        m_model->compute_S_deriv(m_dt);

        // Compute I and R.
        m_model->compute_I_and_R(m_dt);
    }

    std::cout << "Total population at start of simulation is is "
              << m_model->populations.get_value(0)[(Eigen::Index)InfectionState::Susceptible] +
                     m_model->populations.get_value(0)[(Eigen::Index)InfectionState::Infected] +
                     m_model->populations.get_value(0)[(Eigen::Index)InfectionState::Recovered]
              << " and at end "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] +
                     m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] +
                     m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered]
              << std::endl;

    std::cout << "Max number of iterations throughout simulation was " << m_max_number_iterations << std::endl;
}

// In this advance function, we only consider S as in the paper by Messina et al.
void SimulationMessinaExtendedDetailedInit::advance_messina(ScalarType tmax)
{
    // Get index of t0, i.e. index of last time point of given initial values.
    size_t t0_index = m_model->populations.get_num_time_points() - 1;
    std::cout << "t0 index: " << t0_index << std::endl;

    // Set vector with values of transition distribution and parameters, respectively.
    m_model->set_transitiondistribution_vector(m_dt, tmax, t0_index);
    m_model->set_parameter_vectors(m_dt, tmax, t0_index);

    mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {}.", m_model->populations.get_last_time(),
                  tmax, m_dt);

    // Compute S' for t_0,..., t_{n0-1}.
    // Compute S'(0) with forward difference operator due to lack of knowledge of previous values of S.
    // The corresponding flow is then given by -S'.
    m_model->flows.add_time_point(0., TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));
    // m_model->flows.get_value(0)[(Eigen::Index)InfectionTransition::SusceptibleToInfected] =
    //     -(m_model->populations.get_value(1)[(Eigen::Index)InfectionState::Susceptible] -
    //       m_model->populations.get_value(0)[(Eigen::Index)InfectionState::Susceptible]) /
    //     m_dt;

    // Compute S'(t) for t_1,..., t_{n0-1} with backwards difference operator. The corresponding flow is then given by -S'.
    for (size_t i = 1; i < (size_t)m_model->populations.get_num_time_points(); i++) {
        // std::cout << "i: " << i << std::endl;
        m_model->flows.add_time_point(i * m_dt,
                                      TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));
        m_model->flows.get_value(i)[(Eigen::Index)InfectionTransition::SusceptibleToInfected] =
            -(m_model->populations.get_value(i)[(Eigen::Index)InfectionState::Susceptible] -
              m_model->populations.get_value(i - 1)[(Eigen::Index)InfectionState::Susceptible]) /
            m_dt;
    }

    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        // std::cout << "SIR: " << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] << ", "
        //           << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] << ", "
        //           << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] << std::endl;

        // Add new time point to populations.
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));
        m_model->flows.add_time_point(m_model->flows.get_last_time() + m_dt,
                                      Vec::Constant((size_t)InfectionTransition::Count, 0.));

        // Print time.
        if (floating_point_equal(std::remainder(10 * m_model->populations.get_last_time(), tmax), 0., 1e-7)) {
            std::cout << "Time: " << m_model->populations.get_last_time() << std::endl;
        }

        // Compute Susceptibles.
        size_t num_time_points = m_model->populations.get_num_time_points();

        size_t num_iterations = m_model->compute_S(
            m_model->populations.get_value(num_time_points - 2)[(size_t)InfectionState::Susceptible], m_dt, t0_index);

        if (num_iterations > m_max_number_iterations) {
            m_max_number_iterations = num_iterations;
        }

        // Compute S'.
        m_model->compute_S_deriv(m_dt);

        // Compute I and R.
        m_model->compute_I_and_R(m_dt, t0_index);
        // std::cout << "Total population: "
        //           << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] +
        //                  m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] +
        //                  m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered]
        //           << std::endl;
    }

    // std::cout << "SIR: " << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] << ", "
    //           << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] << ", "
    //           << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] << std::endl;

    std::cout << "Total population at start of simulation is "
              << m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Susceptible] +
                     m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Infected] +
                     m_model->populations.get_value(t0_index)[(Eigen::Index)InfectionState::Recovered]
              << " and at end "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] +
                     m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] +
                     m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered]
              << std::endl;

    std::cout << "Max number of iterations throughout simulation was " << m_max_number_iterations << std::endl;
}

// void Simulation::advance(ScalarType tmax)
// {

//     mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {}.", m_model->populations.get_last_time(),
//                   tmax, m_dt);

//     size_t t0_index = m_model->m_finite_difference_order;

//     // Add flows from time 0 until start of the simulation and compute flow from S to I.
//     for (size_t i = 0; i < m_model->populations.get_num_time_points() - t0_index; i++) {
//         m_model->flows.add_time_point(i * m_dt,
//                                       TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));
//         m_model->compute_S_deriv(m_dt, i + t0_index);
//     }

//     std::cout << "Total pop start: " << m_model->populations.get_last_value().sum() << std::endl;
//     while (m_model->populations.get_last_time() < tmax - 1e-10) {

//         // Add new time points to populations and flows.
//         m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
//                                             Vec::Constant((size_t)InfectionState::Count, 0.));
//         m_model->flows.add_time_point(m_model->flows.get_last_time() + m_dt,
//                                       TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));

//         // std::cout << "Current time: " << m_model->populations.get_last_time() << std::endl;

//         // Compute Susceptibles.
//         m_model->compute_S(m_model->populations.get_last_value()[(size_t)InfectionState::Susceptible], m_dt,
//                            m_model->get_totalpop(), t0_index);

//         // Compute flow from S to I via derivative of S.
//         m_model->compute_S_deriv(m_dt);

//         // Compute Infected and Recovered.
//         m_model->compute_I_and_R(m_dt, t0_index, false);
//     }
//     std::cout << "Total pop end: " << m_model->populations.get_last_value().sum() << std::endl;

//     m_model->compute_susceptible_difference(m_dt, t0_index);
// }

// void Simulation::advance2(ScalarType tmax)
// {

//     mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {} with centered difference.",
//                   m_model->populations.get_last_time(), tmax, m_dt);

//     size_t t0_index = m_model->m_finite_difference_order;

//     std::cout << "Total pop start: " << m_model->populations.get_last_value().sum() << std::endl;
//     // Compute only values for S until tmax.
//     while (m_model->populations.get_last_time() < tmax - 1e-10) {

//         // Add new time points to populations and flows.
//         m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
//                                             Vec::Constant((size_t)InfectionState::Count, 0.));

//         // Compute Susceptibles.
//         m_model->compute_S(m_model->populations.get_last_value()[(size_t)InfectionState::Susceptible], m_dt,
//                            m_model->get_totalpop(), t0_index);
//     }

//     // Then compute derivative of S based on centered finite difference and with this I and R.
//     for (size_t i = 0; i < size_t(m_model->populations.get_num_time_points() - t0_index - 2); i++) {
//         m_model->flows.add_time_point(i * m_dt,
//                                       TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));
//         // std::cout << "Time in flows: " << m_model->flows.get_last_time() << std::endl;

//         // Compute flow from S to I via derivative of S.
//         size_t time_point_index_populations = m_model->flows.get_num_time_points() + t0_index - 1;
//         // std::cout << "Time in pop: " << m_model->populations.get_time(time_point_index_populations) << std::endl;
//         m_model->compute_S_deriv_centered(m_dt, time_point_index_populations);

//         if (i >= m_model->get_gregory_order()) {
//             // Compute Infected and Recovered.
//             size_t time_point_index_flows = m_model->flows.get_num_time_points() - 1;
//             m_model->compute_I_and_R_centered(m_dt, t0_index, time_point_index_flows, false);
//         }
//     }
//     std::cout << "Total pop end: " << m_model->populations.get_last_value().sum() << std::endl;

//     m_model->compute_susceptible_difference(m_dt, t0_index);
// }

// TimeSeries<ScalarType> simulate(ScalarType tmax, ScalarType dt, Model const& m_model)
// {
//     Simulation sim(m_model, dt);
//     sim.advance(tmax);
//     return sim.get_result();
// }

} // namespace isir
} // namespace mio
