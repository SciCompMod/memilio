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
#include "ide_sir/simulation_s_groundtruth.h"
#include "memilio/config.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/time_series.h"
#include <cmath>
#include <iostream>

namespace mio
{
namespace isir
{

using Vec = mio::TimeSeries<ScalarType>::Vector;

void SimulationSusceptibleGroundtruth::advance(ScalarType tmax, size_t t0_index)
{
    // // Get index of t0, i.e. index of last time point of given initial values.
    // size_t t0_index = m_model->populations.get_num_time_points() - 1;

    // Set vector with values of transition distribution and parameters, respectively.
    m_model->set_transitiondistribution_vector(m_dt, tmax, t0_index);
    m_model->set_parameter_vectors(m_dt, tmax, t0_index);

    mio::log_info("Simulating I and R from t0 = {} until tmax = {} with dt = {}.",
                  m_model->populations.get_time(t0_index), tmax, m_dt);

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

    // Compute S' as well as I and R.
    // size_t time_point_index = t0_index;
    for (size_t time_point_index = t0_index; time_point_index < (size_t)m_model->flows.get_num_time_points();
         time_point_index++) {

        // Compute I and R.
        m_model->compute_I_and_R(m_dt, time_point_index);
    }

    std::cout << "SIR: " << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] << ", "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered] << std::endl;

    std::cout << "Total population at start of simulation is "
              << m_model->populations.get_value(0)[(Eigen::Index)InfectionState::Susceptible] +
                     m_model->populations.get_value(0)[(Eigen::Index)InfectionState::Infected] +
                     m_model->populations.get_value(0)[(Eigen::Index)InfectionState::Recovered]
              << " and at end "
              << m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] +
                     m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Infected] +
                     m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Recovered]
              << std::endl;
}

} // namespace isir
} // namespace mio
