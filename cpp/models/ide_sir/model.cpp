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
ModelMessinaExtendedDetailedInit::ModelMessinaExtendedDetailedInit(TimeSeries<ScalarType>&& populations_init,
                                                                   ScalarType N_init, size_t gregory_order,
                                                                   size_t finite_difference_order)
    : parameters{Parameters()}
    , populations{std::move(populations_init)}
    , flows{TimeSeries<ScalarType>((size_t)InfectionTransition::Count)}
    , m_N{N_init}
    , m_gregory_order(gregory_order)
    , m_finite_difference_order(finite_difference_order)
{
    assert(m_gregory_order > 0);
    assert(m_finite_difference_order > 0);
}

void ModelMessinaExtendedDetailedInit::set_transitiondistribution_vector(ScalarType dt, ScalarType tmax,
                                                                         size_t t0_index)
{
    m_transitiondistribution_vector = std::vector<ScalarType>(t0_index + size_t(std::ceil(tmax / dt)) + 1, 0.);

    for (size_t i = 0; i <= t0_index + size_t(std::ceil(tmax / dt)); i++) {
        ScalarType state_age = (ScalarType)i * dt;
        m_transitiondistribution_vector[i] =
            parameters.get<TransitionDistributions>()[(size_t)InfectionTransition::InfectedToRecovered].eval(state_age);
    }
}

void ModelMessinaExtendedDetailedInit::set_parameter_vectors(ScalarType dt, ScalarType tmax, size_t t0_index)
{
    m_transmissionproboncontact_vector = std::vector<ScalarType>(t0_index + size_t(std::ceil(tmax / dt)) + 1, 0.);
    m_riskofinffromsymptomatic_vector  = std::vector<ScalarType>(t0_index + size_t(std::ceil(tmax / dt)) + 1, 0.);

    for (size_t i = 0; i <= t0_index + size_t(std::ceil(tmax / dt)); i++) {
        ScalarType state_age                  = (ScalarType)i * dt;
        m_transmissionproboncontact_vector[i] = parameters.get<TransmissionProbabilityOnContact>().eval(state_age);
        m_riskofinffromsymptomatic_vector[i]  = parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age);
    }
}

ScalarType ModelMessinaExtendedDetailedInit::sum_part1_weight(size_t n, size_t j)
{
    assert(n > 0);
    assert(j < n);

    ScalarType gregory_weight;

    if (n < m_gregory_order) {
        std::vector<Eigen::MatrixX<ScalarType>> vec_gregoryweights = get_gregoryweights(n);
        Eigen::MatrixX<ScalarType> gregoryWeights_sigma            = vec_gregoryweights[0];

        gregory_weight = gregoryWeights_sigma(0, j);
    }

    else {

        // Depending on the Gregory order and the current time step, we determine the required row index of gregoryWeights_sigma.
        // This is necessary because we implemented gregoryWeights_sigma in a reduced way.

        size_t row_index;

        // If n <= m_gregory_order + m_gregory_order - 2, then the row_index corresponds to n - gregory_order.
        if (n <= m_gregory_order + m_gregory_order - 2) {
            row_index = n - m_gregory_order;
        }
        // Else, for n > m_gregory_order + m_gregory_order - 2, the entries in gregoryWeights_sigma do not change anymore and the
        // corresponding row_index is given by m_gregory_order - 1.
        else {
            row_index = m_gregory_order - 1;
        }
        // The column index only depends on the current index of the sum j.
        size_t column_index = j;

        std::vector<Eigen::MatrixX<ScalarType>> vec_gregoryweights = get_gregoryweights(m_gregory_order);
        Eigen::MatrixX<ScalarType> gregoryWeights_sigma            = vec_gregoryweights[0];

        gregory_weight = gregoryWeights_sigma(row_index, column_index);
    }

    return gregory_weight;
}

ScalarType ModelMessinaExtendedDetailedInit::sum_part2_weight(size_t n, size_t j)
{
    assert(n > 0);
    assert(j <= n);

    ScalarType gregory_weight;

    if (n < m_gregory_order) {
        std::vector<Eigen::MatrixX<ScalarType>> vec_gregoryweights = get_gregoryweights(n);
        Eigen::MatrixX<ScalarType> gregoryWeights_omega            = vec_gregoryweights[1];

        gregory_weight = gregoryWeights_omega(0);
    }

    else {
        size_t weight_index;

        // Depending on the Gregory order and the current time step, we determine the required weight index of gregoryWeights_omega.
        // This is necessary because we implemented gregoryWeights_omega in a reduced way.
        if (n - j <= m_gregory_order) {
            weight_index = n - j;
        }
        else {
            weight_index = m_gregory_order;
        }

        std::vector<Eigen::MatrixX<ScalarType>> vec_gregoryweights = get_gregoryweights(m_gregory_order);
        Eigen::MatrixX<ScalarType> gregoryWeights_omega            = vec_gregoryweights[1];
        gregory_weight                                             = gregoryWeights_omega(weight_index);
    }

    return gregory_weight;
}

ScalarType ModelMessinaExtendedDetailedInit::fixed_point_function(ScalarType susceptibles, ScalarType dt,
                                                                  size_t t0_index, bool use_complement)
{
    // Get the index of the current time step.
    size_t current_time_index = populations.get_num_time_points() - 1;

    // Compute first part of sum where already known initial values of Susceptibles are used.
    ScalarType sum = 0.;

    // Index determining when we switch from one Gregory sum to the other one.
    // TODO: Explain better what the difference between the two Gregory sums is.
    size_t switch_weights_index = std::min(current_time_index, m_gregory_order);

    if (use_complement) {
        for (size_t j = 0; j <= current_time_index; j++) {
            ScalarType gregory_weight = 0.;
            if (j < switch_weights_index) {
                gregory_weight = sum_part1_weight(current_time_index, j);
            }
            else {
                gregory_weight = sum_part2_weight(current_time_index, j);
            }

            ScalarType relevant_susceptibles;
            // If j<n, we take the Susceptibles at the corresponding index that we have already computed.
            if (j < current_time_index) {
                relevant_susceptibles = populations.get_value(j)[(Eigen::Index)InfectionState::Susceptible];
            }
            // In case of j=n, the number of Susceptibles is not already known and stored in populations but is determined
            // by the input to the fixed point iteration.
            else {
                relevant_susceptibles = susceptibles;
            }

            // For each index, the corresponding summand is computed here.
            sum += gregory_weight *

                   m_transmissionproboncontact_vector[current_time_index - j] *
                   m_riskofinffromsymptomatic_vector[current_time_index - j] *
                   (1. - m_transitiondistribution_vector[current_time_index - j]) *
                   (parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                        SimulationTime<ScalarType>((current_time_index - j + t0_index) * dt))(0, 0) -
                    (parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                         SimulationTime<ScalarType>(current_time_index * dt))(0, 0) /
                     m_N) *
                        relevant_susceptibles);
        }
    }
    else {
        for (size_t j = 0; j <= current_time_index; j++) {

            // // Compute inner sum
            // ScalarType inner_sum = 0.;
            // for (size_t k = 0; k <= j; k++) {
            //     ScalarType gregory_weight_inner_sum = 0.;
            //     if (j < switch_weights_index) {
            //         gregory_weight_inner_sum = sum_part1_weight(current_time_index, j);
            //     }
            //     else {
            //         gregory_weight_inner_sum = sum_part2_weight(current_time_index, j);
            //     }

            //     ScalarType relevant_susceptibles;
            //     // If k<n, we take the Susceptibles at the corresponding index that we have already computed.
            //     if (k < current_time_index) {
            //         relevant_susceptibles = populations.get_value(k)[(Eigen::Index)InfectionState::Susceptible];
            //     }
            //     // In case of j=n, the number of Susceptibles is not already known and stored in populations but is determined
            //     // by the input to the fixed point iteration.
            //     else {
            //         relevant_susceptibles = susceptibles;
            //     }

            //     inner_sum += gregory_weight_inner_sum * m_transmissionproboncontact_vector[j - k] *
            //                  m_riskofinffromsymptomatic_vector[j - k] * m_transitiondistribution_vector[j - k] *
            //                  relevant_susceptibles;
            // }

            ScalarType gregory_weight = 0.;
            if (j < switch_weights_index) {
                gregory_weight = sum_part1_weight(current_time_index, j);
            }
            else {
                gregory_weight = sum_part2_weight(current_time_index, j);
            }

            ScalarType relevant_susceptibles;
            // If j<n, we take the Susceptibles at the corresponding index that we have already computed.
            if (j < current_time_index) {
                relevant_susceptibles = populations.get_value(j)[(Eigen::Index)InfectionState::Susceptible];
            }
            // In case of j=n, the number of Susceptibles is not already known and stored in populations but is determined
            // by the input to the fixed point iteration.
            else {
                relevant_susceptibles = susceptibles;
            }

            // For each index, the corresponding summand is computed here.
            ScalarType current_time = populations.get_last_time();
            sum += gregory_weight * m_transmissionproboncontact_vector[current_time_index - j] *
                   m_riskofinffromsymptomatic_vector[current_time_index - j] *
                   m_transitiondistribution_vector[current_time_index - j] *
                   (parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                        SimulationTime<ScalarType>(current_time - (j + t0_index) * dt))(0, 0) -
                    (parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(
                         SimulationTime<ScalarType>(current_time))(0, 0) /
                     m_N) *
                        relevant_susceptibles);

            // - gregory_weight *
            //     (parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(j * dt)(0, 0) -
            //      parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at((j - 1) * dt)(0, 0)) /
            //     m_N * inner_sum
            // std::cout << "inner sum: " << inner_sum << std::endl;
        }
    }
    return populations.get_value(0)[(Eigen::Index)InfectionState::Susceptible] * std::exp(-dt * sum);
}

size_t ModelMessinaExtendedDetailedInit::compute_S(ScalarType s_init, ScalarType dt, size_t t0_index,
                                                   bool use_complement, ScalarType tol, size_t max_iterations)
{
    size_t iter_counter = 0;
    while (iter_counter < max_iterations) {

        ScalarType s_new = fixed_point_function(s_init, dt, t0_index, use_complement);

        if (std::fabs(s_init - s_new) < tol) {
            break;
        }

        s_init = s_new;
        iter_counter++;
    }

    if (iter_counter == max_iterations) {
        std::cout << "Max number of iterations reached without convergence. Results may not be accurate." << std::endl;
    }

    // Set S in corresponding TimeSeries.
    populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = s_init;

    return iter_counter;
}

void ModelMessinaExtendedDetailedInit::compute_S_deriv(ScalarType dt, size_t time_point_index)
{
    // Linear backwards finite difference scheme, flow from S to I is then given by -S'.
    if (std::min(m_finite_difference_order, time_point_index) == 1) {
        flows[time_point_index][(Eigen::Index)InfectionTransition::SusceptibleToInfected] =
            -(populations[time_point_index][(Eigen::Index)InfectionState::Susceptible] -
              populations[time_point_index - 1][(Eigen::Index)InfectionState::Susceptible]) /
            dt;
    }

    // Compute S' with backwards finite difference scheme of second order, flow from S to I is then given by -S'.
    if (std::min(m_finite_difference_order, time_point_index) == 2) {
        flows[time_point_index][(Eigen::Index)InfectionTransition::SusceptibleToInfected] =
            -(3 * populations[time_point_index][(Eigen::Index)InfectionState::Susceptible] -
              4 * populations[time_point_index - 1][(Eigen::Index)InfectionState::Susceptible] +
              1 * populations[time_point_index - 2][(Eigen::Index)InfectionState::Susceptible]) /
            (2 * dt);
    }

    // Compute S' with backwards finite difference scheme of third order, flow from S to I is then given by -S'.
    if (std::min(m_finite_difference_order, time_point_index) == 3) {
        flows[time_point_index][(Eigen::Index)InfectionTransition::SusceptibleToInfected] =
            -(11 * populations[time_point_index][(Eigen::Index)InfectionState::Susceptible] -
              18 * populations[time_point_index - 1][(Eigen::Index)InfectionState::Susceptible] +
              9 * populations[time_point_index - 2][(Eigen::Index)InfectionState::Susceptible] -
              2 * populations[time_point_index - 3][(Eigen::Index)InfectionState::Susceptible]) /
            (6 * dt);
    }

    // Compute S' with backwards finite difference scheme of fourth order, flow from S to I is then given by -S'.
    if (std::min(m_finite_difference_order, time_point_index) == 4) {
        flows[time_point_index][(Eigen::Index)InfectionTransition::SusceptibleToInfected] =
            -(25 * populations[time_point_index][(Eigen::Index)InfectionState::Susceptible] -
              48 * populations[time_point_index - 1][(Eigen::Index)InfectionState::Susceptible] +
              36 * populations[time_point_index - 2][(Eigen::Index)InfectionState::Susceptible] -
              16 * populations[time_point_index - 3][(Eigen::Index)InfectionState::Susceptible] +
              3 * populations[time_point_index - 4][(Eigen::Index)InfectionState::Susceptible]) /
            (12 * dt);
    }
}

void ModelMessinaExtendedDetailedInit::compute_S_deriv(ScalarType dt)
{
    // Use the number of time points to determine time_point_index, hence we are calculating S deriv for last time point.
    size_t time_point_index = flows.get_num_time_points() - 1;
    compute_S_deriv(dt, time_point_index);
}

void ModelMessinaExtendedDetailedInit::compute_I_and_R(ScalarType dt, size_t time_point_index, bool use_complement)
{
    // Index determining when we switch from one Gregory sum to the other one.
    // TODO: Explain better what the difference between the two Gregory sums is.
    size_t switch_weights_index = std::min(time_point_index, m_gregory_order);

    ScalarType sum_infected  = 0.;
    ScalarType sum_recovered = 0.;

    if (use_complement) {

        // Add first part of sum.
        for (size_t j = 0; j < switch_weights_index; j++) {
            ScalarType gregory_weight = sum_part1_weight(time_point_index, j);

            // For each index, the corresponding summand is computed here.
            sum_infected += gregory_weight * (1. - m_transitiondistribution_vector[time_point_index - j]) *
                            flows.get_value(j)[(Eigen::Index)InfectionTransition::SusceptibleToInfected];
            sum_recovered += gregory_weight * m_transitiondistribution_vector[time_point_index - j] *
                             flows.get_value(j)[(Eigen::Index)InfectionTransition::SusceptibleToInfected];
        }

        // Add second part of sum.
        for (size_t j = switch_weights_index; j <= time_point_index; j++) {
            ScalarType gregory_weight = sum_part2_weight(time_point_index, j);

            // For each index, the corresponding summand is computed here.
            sum_infected += gregory_weight * (1. - m_transitiondistribution_vector[time_point_index - j]) *
                            flows.get_value(j)[(Eigen::Index)InfectionTransition::SusceptibleToInfected];
            sum_recovered += gregory_weight * m_transitiondistribution_vector[time_point_index - j] *
                             flows.get_value(j)[(Eigen::Index)InfectionTransition::SusceptibleToInfected];
        }

        populations[time_point_index][(Eigen::Index)InfectionState::Infected] =
            (1. - m_transitiondistribution_vector[time_point_index]) *
                populations.get_value(0)[(Eigen::Index)InfectionState::Infected] +
            dt * sum_infected;

        populations[time_point_index][(Eigen::Index)InfectionState::Recovered] =
            populations.get_value(0)[(Eigen::Index)InfectionState::Recovered] +
            (m_transitiondistribution_vector[time_point_index]) *
                populations.get_value(0)[(Eigen::Index)InfectionState::Infected] +
            dt * sum_recovered;
    }

    else {
        // Add first part of sum.
        for (size_t j = 0; j < switch_weights_index; j++) {
            ScalarType gregory_weight = sum_part1_weight(time_point_index, j);

            // For each index, the corresponding summand is computed here.
            sum_infected += gregory_weight * m_transitiondistribution_vector[time_point_index - j] *
                            flows.get_value(j)[(Eigen::Index)InfectionTransition::SusceptibleToInfected];
            sum_recovered += gregory_weight * (1. - m_transitiondistribution_vector[time_point_index - j]) *
                             flows.get_value(j)[(Eigen::Index)InfectionTransition::SusceptibleToInfected];
        }

        // Add second part of sum.
        for (size_t j = switch_weights_index; j <= time_point_index; j++) {
            ScalarType gregory_weight = sum_part2_weight(time_point_index, j);

            // For each index, the corresponding summand is computed here.
            sum_infected += gregory_weight * m_transitiondistribution_vector[time_point_index - j] *
                            flows.get_value(j)[(Eigen::Index)InfectionTransition::SusceptibleToInfected];
            sum_recovered += gregory_weight * (1. - m_transitiondistribution_vector[time_point_index - j]) *
                             flows.get_value(j)[(Eigen::Index)InfectionTransition::SusceptibleToInfected];
        }

        populations[time_point_index][(Eigen::Index)InfectionState::Infected] =
            m_transitiondistribution_vector[time_point_index] *
                populations.get_value(0)[(Eigen::Index)InfectionState::Infected] +
            dt * sum_infected;

        // populations[time_point_index][(Eigen::Index)InfectionState::Recovered] =
        //     populations.get_value(0)[(Eigen::Index)InfectionState::Recovered] +
        //     (1. - m_transitiondistribution_vector[time_point_index]) *
        //         populations.get_value(0)[(Eigen::Index)InfectionState::Infected] +
        //     populations[0][(Eigen::Index)InfectionState::Infected] + dt * sum_recovered;
        populations[time_point_index][(Eigen::Index)InfectionState::Recovered] =
            populations.get_value(0)[(Eigen::Index)InfectionState::Recovered] +
            (1. - m_transitiondistribution_vector[time_point_index]) *
                populations.get_value(0)[(Eigen::Index)InfectionState::Infected] +
            populations[0][(Eigen::Index)InfectionState::Susceptible] -
            populations[time_point_index][(Eigen::Index)InfectionState::Susceptible] - dt * sum_infected;
    }
}

void ModelMessinaExtendedDetailedInit::compute_I_and_R(ScalarType dt, bool use_complement)
{
    // Use the number of time points to determine time_point_index, hence we are calculating I and R for last
    // time point of flows.
    size_t time_point_index = flows.get_num_time_points() - 1;
    compute_I_and_R(dt, time_point_index, use_complement);
}

} // namespace isir
} // namespace mio
