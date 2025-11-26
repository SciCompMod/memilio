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
#include "model.h"
#include "gregory_weights.h"
#include "infection_state.h"
#include "memilio/config.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/time_series.h"
#include <Eigen/src/Core/util/Meta.h>
#include <boost/numeric/odeint/util/ublas_wrapper.hpp>
#include <cmath>

namespace mio
{
namespace isir
{
ModelAnalytical::ModelAnalytical(TimeSeries<ScalarType>&& populations_init, size_t gregory_order)
    : populations{std::move(populations_init)}
    , m_gregory_order(gregory_order)
{
    assert(m_gregory_order > 0);
}

ScalarType ModelAnalytical::sum_part1_weight(size_t n, size_t j)
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

ScalarType ModelAnalytical::sum_part2_weight(size_t n, size_t j)
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

ScalarType ModelAnalytical::fixed_point_function(ScalarType susceptibles, ScalarType dt)
{
    // Get the index of the current time step.
    size_t current_time_index = populations.get_num_time_points() - 1;

    // Compute first part of sum where already known initial values of Susceptibles are used.
    ScalarType sum_first_integral = 0.;

    // Index determining when we switch from one Gregory sum to the other one.
    // TODO: Explain better what the difference between the two Gregory sums is.
    size_t switch_weights_index = std::min(current_time_index, m_gregory_order);

    for (size_t j = 0; j < switch_weights_index; j++) {
        ScalarType gregory_weight = sum_part1_weight(current_time_index, j);

        // For each index, the corresponding summand is computed here.
        sum_first_integral += gregory_weight * (current_time_index - j) * dt *
                              populations.get_value(j)[(Eigen::Index)InfectionState::Susceptible];
    }

    // Compute second part of sum where simulated values of Susceptibles are used as well as the given value for the
    // current time step from the fixed point iteration.

    // In this loop, we compute all summands for j=n0,...,n.
    for (size_t j = switch_weights_index; j <= current_time_index; j++) {

        ScalarType gregory_weight = sum_part2_weight(current_time_index, j);

        // For each index, the corresponding summand is computed here.

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

        sum_first_integral += gregory_weight * (current_time_index - j) * dt * relevant_susceptibles;
    }

    return 1 + dt * sum_first_integral;
}

size_t ModelAnalytical::compute_S(ScalarType s_init, ScalarType dt, ScalarType tol, size_t max_iterations)
{
    size_t iter_counter = 0;

    while (iter_counter < max_iterations) {

        // std::cout << "S init: " << s_init << std::endl;
        ScalarType s_new = fixed_point_function(s_init, dt);

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

} // namespace isir
} // namespace mio
