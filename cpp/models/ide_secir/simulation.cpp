/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "ide_secir/parameters.h"
#include "ide_secir/infection_state.h"
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
    mio::log_info("Simulating IDE-SECIR until t={} with dt = {}.", tmax, m_dt);
    m_model->initialize(m_dt);

    // for every time step:
    while (m_model->m_transitions.get_last_time() < tmax - m_dt / 2) {

        m_model->m_transitions.add_time_point(m_model->m_transitions.get_last_time() + m_dt);
        m_model->m_populations.add_time_point(m_model->m_populations.get_last_time() + m_dt);

        // compute_S:
        m_model->compute_susceptibles(m_dt);

        // compute flows:
        m_model->flows_current_timestep(m_dt);

        // compute D
        m_model->compute_deaths();

        // compute m_forceofinfection (only used for calculation of S and sigma_S^E in the next timestep!):
        m_model->update_forceofinfection(m_dt);

        // compute remaining compartments from flows
        m_model->other_compartments_current_timestep(m_dt);
        m_model->compute_recovered();
    }
}

void Simulation::print_transitions() const
{
    // print transitions after simulation
    std::cout << "# time  |  S -> E  |  E - > C  |  C -> I  |  C -> R  |  I -> H  |  I -> R  |  H -> U  |  H -> R  |  "
                 "U -> D  |  U -> R  "
              << std::endl;
    for (Eigen::Index i = 0; i < m_model->m_transitions.get_num_time_points(); ++i) {
        std::cout << m_model->m_transitions.get_time(i);
        for (Eigen::Index j = 0; j < m_model->m_transitions.get_num_elements(); ++j) {
            std::cout << "  |  " << std::fixed << std::setprecision(8) << m_model->m_transitions[i][j];
        }
        std::cout << "\n" << std::endl;
    }
}

void Simulation::print_compartments() const
{
    // print compartments after simulation
    std::cout << "# time  |  S  |  E  |  C  |  I  |  H  |  U  |  R  |  D  |" << std::endl;
    for (Eigen::Index i = 0; i < m_model->m_populations.get_num_time_points(); ++i) {
        std::cout << m_model->m_populations.get_time(i);
        for (Eigen::Index j = 0; j < m_model->m_populations.get_num_elements(); ++j) {
            std::cout << "  |  " << std::fixed << std::setprecision(8) << m_model->m_populations[i][j];
        }
        std::cout << "\n" << std::endl;
    }
}

TimeSeries<ScalarType> simulate(ScalarType t0, ScalarType tmax, ScalarType dt, Model const& m_model)
{
    m_model.check_constraints(dt);
    Simulation sim(m_model, t0, dt);
    sim.advance(tmax);
    return sim.get_result();
}

} // namespace isecir
} // namespace mio
