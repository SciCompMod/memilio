/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Khoa Nguyen
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
#include "abm/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/mioomp.h"
#include "memilio/utils/random_number_generator.h"
#include <random>

namespace mio
{
namespace abm
{

Simulation::Simulation(TimePoint t, World&& world)
    : m_world(std::move(world))
    , m_result(Eigen::Index(InfectionState::Count))
    , m_t(t)
    , m_dt(hours(1))
{
    initialize_locations(m_t);
}

void Simulation::initialize_locations(TimePoint t)
{
    for (auto& location : m_world.get_locations()) {
        location.initialize_subpopulations(t);
    }
}

void Simulation::advance(TimePoint tmax)
{
    //log initial system state
    initialize_locations(m_t);
    store_result_at(m_t);
    while (m_t < tmax) {
        evolve_world(tmax);
        store_result_at(m_t);
    }
}

void Simulation::evolve_world(TimePoint tmax)
{
    auto dt = std::min(m_dt, tmax - m_t);
    m_world.evolve(m_t, dt);
    m_t += m_dt;
}

void Simulation::store_result_at(TimePoint t)
{
    m_result.add_time_point(t.days());
    m_result.get_last_value().setZero();
    PRAGMA_OMP(parallel)
    {
        Eigen::VectorXd sum = Eigen::VectorXd::Zero(m_result.get_num_elements());
        PRAGMA_OMP(for)
        for (auto i = size_t(0); i < m_world.get_locations().size(); ++i) {
            auto&& location = m_world.get_locations()[i];
            sum += location.get_subpopulations().get_last_value().cast<ScalarType>();
        }
        PRAGMA_OMP(critical)
        {
            m_result.get_last_value() += sum;
        }
    }
}

} // namespace abm
} // namespace mio
