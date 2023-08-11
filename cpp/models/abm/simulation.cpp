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
#include "abm/time.h"
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

void Simulation::initialize_locations(TimePoint)
{
    // for (auto& location : m_world.get_locations()) {
    //     location.initialize_subpopulations(t);
    // }
    m_world.prepare();
}

void Simulation::advance(TimePoint tmax)
{
    //log initial system state
    initialize_locations(m_t);
    store_result_at(m_t);
    while (m_t < tmax) {
        auto dt = evolve_world(tmax);
        m_t += dt;
        if (m_t.time_since_midnight() < dt || m_t == tmax) {
            store_result_at(m_t);
        }
    }
}

TimeSpan Simulation::evolve_world(TimePoint tmax)
{
    auto dt = std::min(m_dt, tmax - m_t);
    m_world.evolve(m_t, dt);
    return dt;
}

void Simulation::store_result_at(TimePoint t)
{
    m_result.add_time_point(t.days());
    m_result.get_last_value().setZero();

    //Use a manual parallel reduction to sum up the subpopulations
    //The reduction clause of `omp parallel for` doesn't work well for `Eigen::VectorXd`
    PRAGMA_OMP(parallel)
    {
        //thread local sum of subpopulations, computed in parallel
        //TODO: maybe we can use atomic increments instead of the reduction if necessary for scaling?
        Eigen::VectorXd sum = Eigen::VectorXd::Zero(m_result.get_num_elements());
        auto persons = m_world.get_persons();
        
        PRAGMA_OMP(for schedule(dynamic, 50)) //static?
        for (auto i = size_t(0); i < persons.size(); ++i) {
            auto&& person = persons[i];
            sum[Eigen::Index(person.get_infection_state(t))] += 1;
        }
        
        //synchronized total sum
        PRAGMA_OMP(critical)
        {
            m_result.get_last_value() += sum;
        }
    }
}

} // namespace abm
} // namespace mio
