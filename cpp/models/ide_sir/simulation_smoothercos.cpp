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
#include "ide_sir/model_smoothercos.h"
#include "ide_sir/simulation_smoothercos.h"
#include "memilio/config.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/time_series.h"
#include <bit>
#include <cmath>

namespace mio
{
namespace isir
{

using Vec = mio::TimeSeries<ScalarType>::Vector;

void SimulationSmootherCos::advance(ScalarType tmax)
{

    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        // Add new time point to populations.
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));

        m_model->compute_S();
        m_model->compute_S_deriv_analytical(m_model->populations.get_last_time());
        m_model->compute_S_deriv(m_dt);
    }
}

void SimulationSmootherCos::advance_smoothstep(ScalarType tmax)
{

    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        // Add new time point to populations.
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));

        m_model->compute_S_smoothstep(m_model->populations.get_last_time());
        m_model->compute_smoothstep_deriv_analytical(m_model->populations.get_last_time());
        m_model->compute_smoothstep_deriv_numerical(m_model->populations.get_last_time(), m_dt);
    }
}

} // namespace isir
} // namespace mio
