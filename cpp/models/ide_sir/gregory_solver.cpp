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
#include "ide_sir/gregory_solver.h"
#include "ide_sir/infection_state.h"
#include "ide_sir/parameters.h"
#include "memilio/utils/time_series.h"

namespace mio
{
namespace isir
{

GregorySolver::GregorySolver(TimeSeries<ScalarType> pop, Parameters model_parameters, size_t gregory_order)
    : m_pop(pop)
    , m_model_parameters(model_parameters)
    , m_gregory_order(gregory_order)
{
}

ScalarType GregorySolver::fixed_point_function(ScalarType s, TimeSeries<ScalarType> populations, ScalarType dt,
                                               Parameters model_parameters, ScalarType N)
{
    std::cout << "Time: " << populations.get_last_time() << std::endl;

    Eigen::MatrixX<ScalarType> gregoryWeights_sigma(m_gregory_order, m_gregory_order);
    Eigen::VectorX<ScalarType> gregoryWeights_omega(m_gregory_order + 1);
    size_t scale_gregory_weights;
    if (m_gregory_order == 2) {
        scale_gregory_weights = 12.;
        gregoryWeights_sigma << 5. / scale_gregory_weights, 14. / scale_gregory_weights, 5. / scale_gregory_weights,
            13. / scale_gregory_weights;
        gregoryWeights_omega << 5. / scale_gregory_weights, 13. / scale_gregory_weights, 12. / scale_gregory_weights;
    }
    if (m_gregory_order == 3) {
        scale_gregory_weights = 24.;
        gregoryWeights_sigma << 9. / scale_gregory_weights, 27. / scale_gregory_weights, 27. / scale_gregory_weights,
            9. / scale_gregory_weights, 28. / scale_gregory_weights, 22. / scale_gregory_weights,
            9. / scale_gregory_weights, 28. / scale_gregory_weights, 23. / scale_gregory_weights;
        gregoryWeights_omega << 9. / scale_gregory_weights, 28. / scale_gregory_weights, 23. / scale_gregory_weights,
            24. / scale_gregory_weights;
    }

    ScalarType S_init = populations.get_value(0)[(Eigen::Index)InfectionState::Susceptible];

    ScalarType prefactor = dt * model_parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(0)(0, 0) *
                           model_parameters.get<TransmissionProbabilityOnContact>().eval(0);

    ScalarType sum_init = gregoryWeights_omega(0) * (N - s) / N;
    // std::cout << "gregory weight: " << gregoryWeights_omega( 0) * scale_gregory_weights << std::endl;

    size_t num_time_points = populations.get_num_time_points();

    ScalarType sum_part1 = 0;
    size_t row_index, column_index;
    for (size_t j = 0; j < m_gregory_order; j++) {
        if (num_time_points - m_gregory_order < m_gregory_order) {
            row_index = num_time_points - m_gregory_order;
        }
        else {
            row_index = m_gregory_order - 1;
        }
        column_index = j;

        sum_part1 +=
            gregoryWeights_sigma(row_index, column_index) *
            model_parameters.get<TransitionDistributions>()[(Eigen::Index)InfectionTransition::InfectedToRecovered]
                .eval(num_time_points - j) *
            (N - populations.get_value(j)[(Eigen::Index)InfectionState::Susceptible]);
        // std::cout << "Gregory weight in sum1: " << gregoryWeights_sigma(row_index, column_index) * scale_gregory_weights
        //           << std::endl;
    }

    ScalarType sum_part2 = 0;
    size_t weight_index;
    for (size_t j = m_gregory_order; j < num_time_points; j++) {
        if (num_time_points - j <= m_gregory_order) {
            weight_index = num_time_points - j;
        }
        else {
            weight_index = m_gregory_order;
        }

        sum_part2 +=
            gregoryWeights_omega(weight_index) *
            model_parameters.get<TransitionDistributions>()[(Eigen::Index)InfectionTransition::InfectedToRecovered]
                .eval(num_time_points - j) *
            (N - populations.get_value(j)[(Eigen::Index)InfectionState::Susceptible]);
        // std::cout << "Gregory weight in sum2: " << gregoryWeights_omega(weight_index) * scale_gregory_weights
        //           << std::endl;
    }

    return S_init * std::exp(-prefactor * (sum_init + sum_part1 + sum_part2));
}

void GregorySolver::compute_S(ScalarType s_guess, TimeSeries<ScalarType> populations, ScalarType dt,
                              Parameters model_parameters, ScalarType N, ScalarType tol, size_t max_iterations)
{
    size_t iter_counter = 0;
    while (iter_counter < max_iterations) {

        ScalarType s_estimated = fixed_point_function(s_guess, populations, dt, model_parameters, N);

        if (std::fabs(s_guess - s_estimated) < tol) {
            break;
        }

        s_guess = s_estimated;
        iter_counter++;
    }
    m_pop.get_last_value()[(Eigen::Index)InfectionState::Susceptible] = s_guess;
}

} // namespace isir
} // namespace mio
