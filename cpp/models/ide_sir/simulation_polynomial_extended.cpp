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
#include "ide_sir/model_polynomial_extended.h"
#include "ide_sir/simulation_polynomial_extended.h"
#include "memilio/config.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/time_series.h"
#include <cmath>

namespace mio
{
namespace isir
{

using Vec = mio::TimeSeries<ScalarType>::Vector;

void SimulationPolynomialExtended::advance_single(ScalarType tmax, ScalarType alpha)
{
    // Get index of t0, i.e. index of last time point of given initial values.
    size_t t0_index = m_model->populations.get_num_time_points() - 1;

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
    }

    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        std::cout << "Time pop: " << m_model->populations.get_last_time() << std::endl;

        // Add new time point to populations.
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));

        // Compute Susceptibles.
        m_model->approx_single_integral(m_dt, t0_index, alpha);
    }
}

void SimulationPolynomialExtended::advance_double(ScalarType tmax, ScalarType alpha)
{
    // Get index of t0, i.e. index of last time point of given initial values.
    size_t t0_index = m_model->populations.get_num_time_points() - 1;

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
    }

    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        std::cout << "Time pop: " << m_model->populations.get_last_time() << std::endl;

        // Add new time point to populations.
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));

        // Compute Susceptibles.
        m_model->approx_double_integral(m_dt, t0_index, alpha);
    }
}

} // namespace isir
} // namespace mio
