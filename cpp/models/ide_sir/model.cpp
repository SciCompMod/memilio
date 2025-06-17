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
#include "memilio/config.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/time_series.h"
#include <cmath>

namespace mio
{
namespace isir
{
// Gregory weights corresponding to gregory_order where gregory_order corresponds to n0. The expected order of
// convergence is gregory_order+1.
std::vector<Eigen::MatrixX<ScalarType>> get_gregoryweights(size_t gregory_order)
{
    Eigen::MatrixX<ScalarType> gregoryWeights_sigma(gregory_order, gregory_order);
    Eigen::VectorX<ScalarType> gregoryWeights_omega(gregory_order + 1);
    ScalarType scale_gregory_weights;

    if (gregory_order == 1) {
        scale_gregory_weights = 2.;
        gregoryWeights_sigma << 1 / scale_gregory_weights;
        gregoryWeights_omega << 1. / scale_gregory_weights, 2. / scale_gregory_weights;
    }

    if (gregory_order == 2) {
        scale_gregory_weights = 12.;
        gregoryWeights_sigma << 5. / scale_gregory_weights, 14. / scale_gregory_weights, 5. / scale_gregory_weights,
            13. / scale_gregory_weights;
        gregoryWeights_omega << 5. / scale_gregory_weights, 13. / scale_gregory_weights, 12. / scale_gregory_weights;
    }

    if (gregory_order == 3) {
        scale_gregory_weights = 24.;
        gregoryWeights_sigma << 9. / scale_gregory_weights, 27. / scale_gregory_weights, 27. / scale_gregory_weights,
            9. / scale_gregory_weights, 28. / scale_gregory_weights, 22. / scale_gregory_weights,
            9. / scale_gregory_weights, 28. / scale_gregory_weights, 23. / scale_gregory_weights;
        gregoryWeights_omega << 9. / scale_gregory_weights, 28. / scale_gregory_weights, 23. / scale_gregory_weights,
            24. / scale_gregory_weights;
    }

    return {gregoryWeights_sigma, gregoryWeights_omega};
}

using Vec = mio::TimeSeries<ScalarType>::Vector;

ModelMessina::ModelMessina(TimeSeries<ScalarType>&& populations_init, ScalarType N_init, size_t gregory_order)
    : parameters{Parameters()}
    , populations{std::move(populations_init)}
    , m_N{N_init}
    , m_gregory_order(gregory_order)

{
}

void ModelMessina::set_transitiondistribution_vector(ScalarType dt, ScalarType tmax)
{
    m_transitiondistribution_vector = std::vector<ScalarType>(size_t(std::round(tmax / dt)) + 1, 0.);

    for (size_t i = 0; i <= size_t(std::round(tmax / dt)); i++) {
        ScalarType state_age = (ScalarType)i * dt;
        m_transitiondistribution_vector[i] =
            parameters.get<TransitionDistributions>()[(size_t)InfectionTransition::InfectedToRecovered].eval(state_age);
    }
}

void ModelMessina::set_parameter_vectors(ScalarType dt, ScalarType tmax)
{
    m_transmissionproboncontact_vector = std::vector<ScalarType>(size_t(std::round(tmax / dt)) + 1, 0.);
    m_riskofinffromsymptomatic_vector  = std::vector<ScalarType>(size_t(std::round(tmax / dt)) + 1, 0.);

    for (size_t i = 0; i <= size_t(std::round(tmax / dt)); i++) {
        ScalarType state_age                  = (ScalarType)i * dt;
        m_transmissionproboncontact_vector[i] = parameters.get<TransmissionProbabilityOnContact>().eval(state_age);
        m_riskofinffromsymptomatic_vector[i]  = parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age);
    }
}

ScalarType ModelMessina::sum_part1_term(size_t n, size_t j, ScalarType input)
{
    // Depending on the Gregory order and the current time step, we determine the required row index of gregoryWeights_sigma.
    // This is necessary because we implemented gregoryWeights_sigma in a reduced way.

    size_t row_index, column_index;

    // If n == gregory_order, then the row_index corresponds to n - gregory_order.
    if (n <= m_gregory_order + m_gregory_order - 2) {
        row_index = n - m_gregory_order;
    }
    // Else, for n >= m_gregory_order - 1, the entries in gregoryWeights_sigma do not change anymore and the
    // corresponding row_index is given by m_gregory_order - 1.
    else {
        row_index = m_gregory_order - 1;
    }
    // The column index only depends on the current index of the sum j.
    column_index = j;

    std::vector<Eigen::MatrixX<ScalarType>> vec_gregoryweights = get_gregoryweights(m_gregory_order);
    Eigen::MatrixX<ScalarType> gregoryWeights_sigma            = vec_gregoryweights[0];

    ScalarType sum_part1_term;

    sum_part1_term = gregoryWeights_sigma(row_index, column_index) * m_transmissionproboncontact_vector[n - j] *
                     m_riskofinffromsymptomatic_vector[n - j] * m_transitiondistribution_vector[n - j] * input;

    return sum_part1_term;
}

ScalarType ModelMessina::sum_part2_term(size_t n, size_t j, ScalarType input)
{
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

    ScalarType sum_part2_term;

    sum_part2_term = gregoryWeights_omega(weight_index) * m_transmissionproboncontact_vector[n - j] *
                     m_riskofinffromsymptomatic_vector[n - j] * m_transitiondistribution_vector[n - j] * input;

    return sum_part2_term;
}

ScalarType ModelMessina::fixed_point_function(ScalarType susceptibles, ScalarType dt)
{
    // Get Gregory weights.
    // TODO: Do this in a private function of Simulation class? Or solver class?
    std::vector<Eigen::MatrixX<ScalarType>> vec_gregoryweights = get_gregoryweights(m_gregory_order);
    Eigen::MatrixX<ScalarType> gregoryWeights_sigma            = vec_gregoryweights[0];
    Eigen::MatrixX<ScalarType> gregoryWeights_omega            = vec_gregoryweights[1];

    // Get S0.
    ScalarType S0 = populations.get_value(0)[(Eigen::Index)InfectionState::Susceptible];

    // Get current time.
    ScalarType current_time = populations.get_last_time();

    // Get the index of the current time step.
    size_t n = populations.get_num_time_points() - 1;

    // Compute first part of sum where already known initial values of Susceptibles are used.
    ScalarType sum_part1 = 0;

    for (size_t j = 0; j < m_gregory_order; j++) {

        // For each index, the corresponding summand is computed here.
        sum_part1 += sum_part1_term(n, j, (m_N - populations.get_value(j)[(Eigen::Index)InfectionState::Susceptible]));
    }
    // std::cout << "Sum part 1: " << sum_part1 << std::endl;

    // Compute second part of sum where simulated values of Susceptibles are used.
    ScalarType sum_part2 = 0;

    // In this loop, we compute all summands for j=n0,...,n-1. The summand for j=n is added separately below.
    for (size_t j = m_gregory_order; j <= n; j++) {

        // For each index, the corresponding summand is computed here.
        if (j < n) {
            sum_part2 +=
                sum_part2_term(n, j, m_N - populations.get_value(j)[(Eigen::Index)InfectionState::Susceptible]);
        }
        // In case of j=n, the number of Susceptibles is not already known and stored in populations but is determined
        // by the fixed point iteration.
        else {
            sum_part2 += sum_part2_term(n, j, m_N - susceptibles);
        }
    }
    // std::cout << "Sum part 2: " << sum_part2 << std::endl;

    ScalarType prefactor = dt * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(current_time)(0, 0);
    return S0 * std::exp(-prefactor * (sum_part1 + sum_part2));
}

void ModelMessina::compute_S(ScalarType s_init, ScalarType dt, ScalarType tol, size_t max_iterations)
{
    size_t iter_counter = 0;
    while (iter_counter < max_iterations) {

        ScalarType s_estimated = fixed_point_function(s_init, dt);
        // std::cout << "s_estimated: " << s_estimated << std::endl;

        // std::cout << "Diff: " << std::fabs(s_init - s_estimated) << "; tol: " << tol << std::endl;
        if (std::fabs(s_init - s_estimated) < tol) {
            break;
        }

        s_init = s_estimated;
        iter_counter++;
    }

    if (iter_counter == max_iterations) {
        std::cout << "Max number of iterations reached without convergence. Results may not be accurate." << std::endl;
    }

    // std::cout << "Num iterations: " << iter_counter << std::endl;

    // Set S in corresponding TimeSeries.
    populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = s_init;
}

/*********************************************************************************************************************/

Model::Model(TimeSeries<ScalarType>&& populations_init, ScalarType N_init, size_t gregory_order,
             size_t finite_difference_order)
    : parameters{Parameters()}
    , populations{std::move(populations_init)}
    , flows(TimeSeries<ScalarType>((size_t)InfectionTransition::Count))
    , m_finite_difference_order(finite_difference_order)
    , m_N{N_init}
    , m_gregory_order(gregory_order)

{
}

ScalarType Model::get_totalpop() const
{
    return m_N;
}

ScalarType Model::sum_part1_term(size_t row_index, size_t column_index, ScalarType state_age, ScalarType input,
                                 bool recovered)
{

    std::vector<Eigen::MatrixX<ScalarType>> vec_gregoryweights = get_gregoryweights(m_gregory_order);
    Eigen::MatrixX<ScalarType> gregoryWeights_sigma            = vec_gregoryweights[0];

    ScalarType sum_part1_term;
    if (!recovered) {
        sum_part1_term =
            gregoryWeights_sigma(row_index, column_index) *
            parameters.get<TransmissionProbabilityOnContact>().eval(state_age) *
            parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age) *
            parameters.get<TransitionDistributions>()[(Eigen::Index)InfectionTransition::InfectedToRecovered].eval(
                state_age) *
            input;
    }
    else {
        sum_part1_term =
            gregoryWeights_sigma(row_index, column_index) *
            parameters.get<TransmissionProbabilityOnContact>().eval(state_age) *
            parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age) *
            (1 - parameters.get<TransitionDistributions>()[(Eigen::Index)InfectionTransition::InfectedToRecovered].eval(
                     state_age)) *
            input;
    }

    // std::cout << "Gregory weight in sum1: " << gregoryWeights_sigma(row_index, column_index) << std::endl;

    return sum_part1_term;
}

ScalarType Model::sum_part2_term(size_t weight_index, ScalarType state_age, ScalarType input, bool recovered)
{
    std::vector<Eigen::MatrixX<ScalarType>> vec_gregoryweights = get_gregoryweights(m_gregory_order);
    Eigen::MatrixX<ScalarType> gregoryWeights_omega            = vec_gregoryweights[1];

    ScalarType sum_part2_term;
    if (!recovered) {
        sum_part2_term =
            gregoryWeights_omega(weight_index) * parameters.get<TransmissionProbabilityOnContact>().eval(state_age) *
            parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age) *
            parameters.get<TransitionDistributions>()[(Eigen::Index)InfectionTransition::InfectedToRecovered].eval(
                state_age) *
            input;
    }
    else {
        sum_part2_term =
            gregoryWeights_omega(weight_index) * parameters.get<TransmissionProbabilityOnContact>().eval(state_age) *
            parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age) *
            (1 - parameters.get<TransitionDistributions>()[(Eigen::Index)InfectionTransition::InfectedToRecovered].eval(
                     state_age)) *
            input;
    }
    // std::cout << "Gregory weight in sum2: " << gregoryWeights_omega(weight_index) << std::endl;

    return sum_part2_term;
}

ScalarType Model::fixed_point_function(ScalarType s, ScalarType dt, ScalarType N, size_t t0_index)
{
    // TODO: Do this in a private function of Simulation class? Or solver class?
    std::vector<Eigen::MatrixX<ScalarType>> vec_gregoryweights = get_gregoryweights(m_gregory_order);
    Eigen::MatrixX<ScalarType> gregoryWeights_sigma            = vec_gregoryweights[0];
    Eigen::MatrixX<ScalarType> gregoryWeights_omega            = vec_gregoryweights[1];

    ScalarType S_init = populations.get_value(t0_index)[(Eigen::Index)InfectionState::Susceptible];

    ScalarType current_time = populations.get_last_time() - dt;

    ScalarType prefactor = dt * parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(current_time)(0, 0);

    ScalarType sum_init = gregoryWeights_omega(0) * parameters.get<TransmissionProbabilityOnContact>().eval(0.) *
                          parameters.get<RiskOfInfectionFromSymptomatic>().eval(0.) * (N - s) / N;
    // std::cout << "Gregory weight in sum_init: " << gregoryWeights_omega(0) << std::endl;

    size_t num_time_points_simulated = populations.get_num_time_points() - t0_index - 2;

    ScalarType sum_part1 = 0;
    size_t row_index, column_index;
    for (size_t j = 0; j < m_gregory_order; j++) {
        if (num_time_points_simulated - m_gregory_order < m_gregory_order) {
            row_index = num_time_points_simulated - m_gregory_order;
        }
        else {
            row_index = m_gregory_order - 1;
        }
        column_index = j;

        ScalarType state_age = (num_time_points_simulated - j) * dt;

        sum_part1 +=
            sum_part1_term(row_index, column_index, state_age,
                           (N - populations.get_value(j + t0_index)[(Eigen::Index)InfectionState::Susceptible]));
    }

    ScalarType sum_part2 = 0;
    size_t weight_index;
    for (size_t j = m_gregory_order; j < num_time_points_simulated; j++) {
        if (num_time_points_simulated - j + t0_index <= m_gregory_order) {
            weight_index = num_time_points_simulated - j + t0_index;
        }
        else {
            weight_index = m_gregory_order;
        }

        ScalarType state_age = (num_time_points_simulated - j) * dt;

        sum_part2 +=
            sum_part2_term(weight_index, state_age,
                           (N - populations.get_value(j + t0_index)[(Eigen::Index)InfectionState::Susceptible]));
    }

    return S_init * std::exp(-prefactor * (sum_init + sum_part1 + sum_part2));
}

void Model::compute_S(ScalarType s_init, ScalarType dt, ScalarType N, size_t t0_index, ScalarType tol,
                      size_t max_iterations)
{
    size_t iter_counter = 0;
    while (iter_counter < max_iterations) {

        ScalarType s_estimated = fixed_point_function(s_init, dt, N, t0_index);

        if (std::fabs(s_init - s_estimated) < tol) {
            break;
        }

        s_init = s_estimated;
        iter_counter++;
    }

    if (iter_counter == max_iterations) {
        std::cout << "Max number of iterations reached without convergence. Results may not be accurate." << std::endl;
    }

    // std::cout << "Num iterations: " << iter_counter << std::endl;

    // Set S.

    populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = s_init;
}

void Model::compute_S_deriv(ScalarType dt, size_t time_point_index)
{
    // Compute S' with finite difference scheme of fourth order, flow from S to I is then given by -S'.
    if (m_finite_difference_order == 4) {
        flows.get_last_value()[(Eigen::Index)InfectionTransition::SusceptibleToInfected] =
            -(3 * populations[time_point_index - 4][(Eigen::Index)InfectionState::Susceptible] -
              16 * populations[time_point_index - 3][(Eigen::Index)InfectionState::Susceptible] +
              36 * populations[time_point_index - 2][(Eigen::Index)InfectionState::Susceptible] -
              48 * populations[time_point_index - 1][(Eigen::Index)InfectionState::Susceptible] +
              25 * populations[time_point_index][(Eigen::Index)InfectionState::Susceptible]) /
            (12 * dt);
    }

    // Linear finite difference scheme
    if (m_finite_difference_order == 1) {
        flows.get_last_value()[(Eigen::Index)InfectionTransition::SusceptibleToInfected] =
            -(populations[time_point_index][(Eigen::Index)InfectionState::Susceptible] -
              populations[time_point_index - 1][(Eigen::Index)InfectionState::Susceptible]) /
            dt;
    }
}

void Model::compute_S_deriv(ScalarType dt)
{
    // Use the number of time points to determine time_point_index, hence we are calculating S deriv for last time point.
    size_t time_point_index = populations.get_num_time_points() - 1;
    compute_S_deriv(dt, time_point_index);
}

void Model::compute_I_and_R(ScalarType dt, size_t t0_index, bool enforce_mass_conservation)
{
    size_t num_time_points_simulated = populations.get_num_time_points() - t0_index;

    // Compute I and R.
    ScalarType sum_part1_I = 0;
    ScalarType sum_part1_R = 0;
    size_t row_index, column_index;
    for (size_t j = 0; j < m_gregory_order; j++) {
        if (num_time_points_simulated - m_gregory_order + t0_index < m_gregory_order) {
            row_index = num_time_points_simulated - m_gregory_order + t0_index;
        }
        else {
            row_index = m_gregory_order - 1;
        }
        column_index = j;

        ScalarType state_age = (num_time_points_simulated - j) * dt;

        sum_part1_I += sum_part1_term(
            row_index, column_index, state_age,
            flows[j][(Eigen::Index)InfectionTransition::SusceptibleToInfected] +
                populations[populations.get_num_time_points() - 2][(Eigen::Index)InfectionState::Infected]);

        sum_part1_R += sum_part1_term(
            row_index, column_index, state_age,
            flows[j][(Eigen::Index)InfectionTransition::SusceptibleToInfected] +
                populations[populations.get_num_time_points() - 2][(Eigen::Index)InfectionState::Infected],
            true);
    }

    ScalarType sum_part2_I = 0;
    ScalarType sum_part2_R = 0;

    size_t weight_index;
    for (size_t j = m_gregory_order; j < num_time_points_simulated; j++) {
        if (num_time_points_simulated - j + t0_index <= m_gregory_order) {
            weight_index = num_time_points_simulated - j + t0_index;
        }
        else {
            weight_index = m_gregory_order;
        }

        ScalarType state_age = (num_time_points_simulated - j) * dt;

        sum_part2_I +=
            sum_part2_term(weight_index, state_age, flows[j][(Eigen::Index)InfectionTransition::SusceptibleToInfected]);

        sum_part2_R += sum_part2_term(weight_index, state_age,
                                      flows[j][(Eigen::Index)InfectionTransition::SusceptibleToInfected], true);
    }

    populations.get_last_value()[(Eigen::Index)InfectionState::Infected] = dt * (sum_part1_I + sum_part2_I);

    populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] =
        populations[t0_index][(Eigen::Index)InfectionState::Recovered] + dt * (sum_part1_R + sum_part2_R);

    if (enforce_mass_conservation) {
        ScalarType mass_conservation_scaling =
            (m_N - populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible]) /
            (populations.get_last_value()[(Eigen::Index)InfectionState::Infected] +
             populations.get_last_value()[(Eigen::Index)InfectionState::Recovered]);
        populations.get_last_value()[(Eigen::Index)InfectionState::Infected] *= mass_conservation_scaling;
        populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] *= mass_conservation_scaling;
    }
}

void Model::compute_S_deriv_centered(ScalarType dt, size_t time_point_index)
{
    // std::cout << "time in populations: " << populations.get_time(time_point_index) << std::endl;
    // Compute S' with a centered finite difference scheme of fourth order, flow from S to I is then given by -S'.
    flows.get_last_value()[(Eigen::Index)InfectionTransition::SusceptibleToInfected] =
        -(1 * populations[time_point_index - 2][(Eigen::Index)InfectionState::Susceptible] -
          8 * populations[time_point_index - 1][(Eigen::Index)InfectionState::Susceptible] +
          8 * populations[time_point_index + 1][(Eigen::Index)InfectionState::Susceptible] -
          1 * populations[time_point_index + 2][(Eigen::Index)InfectionState::Susceptible]) /
        (12 * dt);
}

void Model::compute_I_and_R_centered(ScalarType dt, size_t t0_index, size_t time_point_index,
                                     bool enforce_mass_conservation)
{
    size_t num_time_points_simulated = time_point_index + m_gregory_order;

    // Compute I and R.
    ScalarType sum_part1_I = 0;
    ScalarType sum_part1_R = 0;
    size_t row_index, column_index;
    for (size_t j = 0; j < std::min(m_gregory_order, time_point_index); j++) {
        if (num_time_points_simulated - m_gregory_order < m_gregory_order) {
            row_index = num_time_points_simulated - m_gregory_order;
        }
        else {
            row_index = m_gregory_order - 1;
        }
        column_index = j;

        ScalarType state_age = (num_time_points_simulated - j) * dt;

        sum_part1_I += sum_part1_term(row_index, column_index, state_age,
                                      flows[j][(Eigen::Index)InfectionTransition::SusceptibleToInfected]);

        sum_part1_R += sum_part1_term(row_index, column_index, state_age,
                                      flows[j][(Eigen::Index)InfectionTransition::SusceptibleToInfected], true);
    }

    ScalarType sum_part2_I = 0;
    ScalarType sum_part2_R = 0;

    size_t weight_index;
    for (size_t j = m_gregory_order; j < time_point_index; j++) {
        if (num_time_points_simulated - j <= m_gregory_order) {
            weight_index = num_time_points_simulated - j;
        }
        else {
            weight_index = m_gregory_order;
        }

        ScalarType state_age = (num_time_points_simulated - j) * dt;

        sum_part2_I +=
            sum_part2_term(weight_index, state_age, flows[j][(Eigen::Index)InfectionTransition::SusceptibleToInfected]);

        sum_part2_R += sum_part2_term(weight_index, state_age,
                                      flows[j][(Eigen::Index)InfectionTransition::SusceptibleToInfected], true);
    }

    populations[time_point_index + t0_index][(Eigen::Index)InfectionState::Infected] = dt * (sum_part1_I + sum_part2_I);

    populations[time_point_index + t0_index][(Eigen::Index)InfectionState::Recovered] =
        populations[t0_index][(Eigen::Index)InfectionState::Recovered] + dt * (sum_part1_R + sum_part2_R);

    if (enforce_mass_conservation) {
        ScalarType mass_conservation_scaling =
            (m_N - populations[time_point_index + t0_index][(Eigen::Index)InfectionState::Susceptible]) /
            (populations[time_point_index + t0_index][(Eigen::Index)InfectionState::Infected] +
             populations[time_point_index + t0_index][(Eigen::Index)InfectionState::Recovered]);

        // std::cout << "mass_conservation_scaling: " << mass_conservation_scaling << std::endl;
        populations[time_point_index + t0_index][(Eigen::Index)InfectionState::Infected] *= mass_conservation_scaling;
        populations[time_point_index + t0_index][(Eigen::Index)InfectionState::Recovered] *= mass_conservation_scaling;
    }
}

void Model::compute_susceptible_difference(ScalarType dt, size_t t0_index)
{
    for (Eigen::Index i = 1; i < populations.get_num_time_points(); i++) {
        susceptibles_difference.add_time_point(
            ((ScalarType)i - t0_index) * dt,
            Vec::Constant(1, populations[i][(Eigen::Index)InfectionState::Susceptible] -
                                 populations[i - 1][(Eigen::Index)InfectionState::Susceptible]));
    }
}

} // namespace isir
} // namespace mio
