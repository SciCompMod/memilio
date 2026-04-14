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
#include "ide_sir/model_polynomial_extended.h"
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
ModelPolynomialExtended::ModelPolynomialExtended(TimeSeries<ScalarType>&& populations_init, size_t gregory_order,
                                                 ScalarType N)
    : populations{std::move(populations_init)}
    , flows{TimeSeries<ScalarType>((size_t)InfectionTransition::Count)}
    , m_gregory_order(gregory_order)
    , m_N(N)
{
    assert(m_gregory_order > 0);
    std::cout << "total pop: " << m_N << std::endl;
}

ScalarType ModelPolynomialExtended::sum_part1_weight(size_t n, size_t j)
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

ScalarType ModelPolynomialExtended::sum_part2_weight(size_t n, size_t j)
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

void ModelPolynomialExtended::approx_single_integral(ScalarType dt, size_t t0_index, ScalarType alpha)
{
    unused(t0_index);
    // Get the index of the current time step.
    int current_time_index = populations.get_num_time_points() - 1;
    std::cout << "current time index: " << current_time_index << std::endl;
    // std::cout << "t0 index: " << t0_index << std::endl;

    // Compute first part of sum where already known initial values of Susceptibles are used.
    ScalarType sum = 0.;

    for (int j = 0; j <= (int)current_time_index; j++) {

        ScalarType gregory_weight = 0.;
        int switch_weights_index  = std::min(current_time_index, (int)m_gregory_order);
        if (j < switch_weights_index) {
            gregory_weight = sum_part1_weight(current_time_index, j);
        }
        else {
            gregory_weight = sum_part2_weight(current_time_index, j);
        }

        // For each index, the corresponding summand is computed here.
        sum += dt * gregory_weight *
               ((1 - alpha * current_time_index * dt) * std::pow(j * dt, 4) -
                (1 - alpha * (ScalarType(current_time_index) - ScalarType(j) - ScalarType(t0_index)) * dt) * m_N +
                alpha * 1. / 5. * std::pow(j * dt, 5));
    }

    populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = sum;
}

void ModelPolynomialExtended::approx_double_integral(ScalarType dt, size_t t0_index, ScalarType alpha)
{
    unused(t0_index);
    // Get the index of the current time step.
    int current_time_index = populations.get_num_time_points() - 1;
    std::cout << "current time index: " << current_time_index << std::endl;

    // Compute first part of sum where already known initial values of Susceptibles are used.
    ScalarType sum = 0.;

    for (int j = 0; j <= (int)current_time_index; j++) {
        std::cout << "j: " << j << std::endl;
        // Compute inner sum
        ScalarType inner_sum = 0.;
        if (j > 0) {

            for (int k = 0; k <= j; k++) {
                ScalarType gregory_weight_inner_sum = 0.;
                int switch_weights_index_inner_sum  = std::min(j, (int)m_gregory_order);
                if (k < switch_weights_index_inner_sum) {
                    gregory_weight_inner_sum = sum_part1_weight(j, k);
                }
                else {
                    gregory_weight_inner_sum = sum_part2_weight(j, k);
                }

                inner_sum += gregory_weight_inner_sum * std::pow(k * dt, 4);
            }
        }

        ScalarType gregory_weight = 0.;
        int switch_weights_index  = std::min(current_time_index, (int)m_gregory_order);
        if (j < switch_weights_index) {
            gregory_weight = sum_part1_weight(current_time_index, j);
        }
        else {
            gregory_weight = sum_part2_weight(current_time_index, j);
        }

        // For each index, the corresponding summand is computed here.

        sum += dt * gregory_weight *
               ((1 - alpha * current_time_index * dt) * std::pow(j * dt, 4) -
                (1 - alpha * (ScalarType(current_time_index) - ScalarType(j) - ScalarType(t0_index)) * dt) * m_N +
                dt * alpha * inner_sum);
    }

    populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = sum;
}

} // namespace isir
} // namespace mio
