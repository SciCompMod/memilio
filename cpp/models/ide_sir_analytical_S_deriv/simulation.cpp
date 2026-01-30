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

#include "ide_sir_analytical_S_deriv/infection_state.h"
#include "ide_sir_analytical_S_deriv/model.h"
#include "ide_sir_analytical_S_deriv/simulation.h"
#include "memilio/config.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/time_series.h"
#include <cmath>
#include <cstddef>

namespace mio
{
namespace isir
{

using Vec = mio::TimeSeries<ScalarType>::Vector;

void SimulationAnalyticalSDeriv::advance(ScalarType tmax)
{
    mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {}.", m_model->populations.get_last_time(),
                  tmax, m_dt);

    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        // Print time.
        if (floating_point_equal(std::remainder(10 * m_model->populations.get_last_time(), tmax), 0., 1e-7)) {
            std::cout << "Time pop: " << m_model->populations.get_last_time() << std::endl;
        }

        // Add new time point to populations.
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));

        // Compute Susceptibles.
        size_t num_time_points = m_model->populations.get_num_time_points();

        size_t num_iterations = m_model->compute_S(
            m_model->populations.get_value(num_time_points - 2)[(size_t)InfectionState::Susceptible], m_dt);

        if (num_iterations > m_max_number_iterations) {
            m_max_number_iterations = num_iterations;
        }
    }

    // Compute S' as well as I and R.
    for (size_t i = 0; i < (size_t)m_model->populations.get_num_time_points(); i++) {

        m_model->flows.add_time_point(i * m_dt, Vec::Constant((size_t)InfectionTransition::Count, 0.));

        if (floating_point_equal(std::remainder(10 * m_model->flows.get_last_time(), tmax), 0., 1e-7)) {
            std::cout << "Time flows: " << m_model->flows.get_last_time() << std::endl;
        }

        m_model->compute_S_deriv(m_dt);
    }

    // for (size_t i = 0; i < m_model->get_finite_difference_order() + 1; i++) {
    //     m_model->populations.remove_time_point(i);
    //     m_model->flows.remove_time_point(i);
    // }

    for (size_t i = 0; i < (size_t)m_model->populations.get_num_time_points(); i++) {
        // Compute I and R.
        m_model->compute_I_and_R(m_dt, i);
    }

    std::cout << "SIR: " << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] << std::endl;

    std::cout << "Max number of iterations throughout simulation was " << m_max_number_iterations << std::endl;

    std::cout << std::endl;
}

} // namespace isir
} // namespace mio
