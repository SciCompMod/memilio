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

void Simulation::advance(ScalarType tmax)
{
    mio::log_info("Simulating IDE-SIR from t0 = {} until tmax = {} with dt = {}.", m_model->populations.get_last_time(),
                  tmax, m_dt);

    // For every time step:
    using Vec = mio::TimeSeries<ScalarType>::Vector;
    while (m_model->populations.get_last_time() < tmax - 1e-10) {

        m_solver.compute_S(m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible],
                           m_model->populations, m_dt, m_model->parameters, m_model->get_totalpop());

        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt,
                                            Vec::Constant((size_t)InfectionState::Count, 0.));
        m_model->populations.get_last_value()[(Eigen::Index)InfectionState::Susceptible] =
            m_solver.m_pop.get_last_value()[(Eigen::Index)InfectionState::Susceptible];
        m_solver.m_pop = m_model->populations;
    }
}

TimeSeries<ScalarType> simulate(ScalarType tmax, ScalarType dt, Model const& m_model, size_t gregory_order)
{
    Simulation sim(m_model, dt, gregory_order);
    sim.advance(tmax);
    return sim.get_result();
}

} // namespace isir
} // namespace mio
