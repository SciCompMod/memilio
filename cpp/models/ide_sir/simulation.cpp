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
#include "memilio/utils/time_series.h"

namespace mio
{
namespace isir
{

using Vec = mio::TimeSeries<ScalarType>::Vector;

void Simulation::advance(ScalarType tmax)
{

    mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {}.", m_model->populations.get_last_time(),
                  tmax, m_dt);

    size_t t0_index = m_model->m_finite_difference_order;

    // Add flows from time 0 until start of the simulation and compute flow from S to I.
    for (size_t i = 0; i < m_model->populations.get_num_time_points() - t0_index; i++) {
        m_model->flows.add_time_point(i * m_dt,
                                      TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));
        m_model->compute_S_deriv(m_dt, i + t0_index);
    }

    std::cout << "Total pop: " << m_model->populations.get_last_value().sum() << std::endl;
    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        // Add new time points to populations and flows.
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));
        m_model->flows.add_time_point(m_model->flows.get_last_time() + m_dt,
                                      TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));

        // Compute Susceptibles.
        m_model->compute_S(m_model->populations.get_last_value()[(size_t)InfectionState::Susceptible], m_dt,
                           m_model->get_totalpop(), t0_index);

        // Compute flow from S to I via derivative of S.
        m_model->compute_S_deriv(m_dt);

        // Compute Infected and Recovered.
        m_model->compute_I_and_R(m_dt, t0_index);

        std::cout << "Total pop: " << m_model->populations.get_last_value().sum() << std::endl;
    }
}

void Simulation::advance2(ScalarType tmax)
{

    mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {} with centered difference.",
                  m_model->populations.get_last_time(), tmax, m_dt);

    size_t t0_index = m_model->m_finite_difference_order;

    std::cout << "Total pop: " << m_model->populations.get_last_value().sum() << std::endl;
    // Compute only values for S until tmax.
    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        // Add new time points to populations and flows.
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));

        // Compute Susceptibles.
        m_model->compute_S(m_model->populations.get_last_value()[(size_t)InfectionState::Susceptible], m_dt,
                           m_model->get_totalpop(), t0_index);
    }

    // Then compute derivative of S based on centered finite difference and with this I and R.
    for (size_t i = 0; i < size_t(m_model->populations.get_num_time_points() - t0_index - 2); i++) {
        m_model->flows.add_time_point(i * m_dt,
                                      TimeSeries<ScalarType>::Vector::Constant((size_t)InfectionTransition::Count, 0.));
        // std::cout << "Time in flows: " << m_model->flows.get_last_time() << std::endl;

        // Compute flow from S to I via derivative of S.
        size_t time_point_index_populations = m_model->flows.get_num_time_points() + t0_index - 1;
        // std::cout << "Time in pop: " << m_model->populations.get_time(time_point_index_populations) << std::endl;
        m_model->compute_S_deriv_centered(m_dt, time_point_index_populations);

        if (i >= m_model->get_gregory_order()) {
            // Compute Infected and Recovered.
            size_t time_point_index_flows = m_model->flows.get_num_time_points() - 1;
            m_model->compute_I_and_R_centered(m_dt, t0_index, time_point_index_flows);
        }
        // std::cout << "Total pop: " << m_model->populations[time_point_index_populations].sum() << std::endl;
    }
}

TimeSeries<ScalarType> simulate(ScalarType tmax, ScalarType dt, Model const& m_model)
{
    Simulation sim(m_model, dt);
    sim.advance(tmax);
    return sim.get_result();
}

} // namespace isir
} // namespace mio
