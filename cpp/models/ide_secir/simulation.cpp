/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Martin J Kuehn, Anna Wendler, Lena Ploetzke
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
#include "ide_secir/simulation.h"
#include "ide_secir/model.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include <memory>

namespace mio
{
namespace isecir
{

void Simulation::advance(ScalarType tmax)
{
    mio::log_info("Simulating IDE-SECIR from t0 = {} until tmax = {} with dt = {}.",
                  m_model->transitions.get_last_time(), tmax, m_dt);
    m_model->set_transitiondistributions_support_max(m_dt);
    m_model->set_transitiondistributions_derivative(m_dt);
    m_model->set_transitiondistributions_in_forceofinfection(m_dt);
    m_model->initial_compute_compartments(m_dt);

    // For every time step:
    while (m_model->transitions.get_last_time() < tmax - m_dt / 2) {

        m_model->transitions.add_time_point(m_model->transitions.get_last_time() + m_dt);
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt);

        // compute Susceptibles:
        m_model->compute_susceptibles(m_dt);

        // Compute flows:
        m_model->flows_current_timestep(m_dt);

        // Update remaining compartments:
        m_model->update_compartments();

        // Compute m_forceofinfection (only used for calculation of Susceptibles and flow SusceptibleToExposed in the next timestep!):
        m_model->compute_forceofinfection(m_dt);
    }
}

TimeSeries<ScalarType> simulate(ScalarType tmax, ScalarType dt, Model const& m_model)
{
    m_model.check_constraints(dt);
    Simulation sim(m_model, dt);
    sim.advance(tmax);
    return sim.get_result();
}

} // namespace isecir
} // namespace mio
