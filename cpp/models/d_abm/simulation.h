/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: René Schmieding, Julia Bicker
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

#ifndef MIO_D_ABM_SIMULATION_H
#define MIO_D_ABM_SIMULATION_H

#include "d_abm/model.h"
#include "memilio/config.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/random_number_generator.h"

namespace mio
{

namespace dabm
{

/**
 * @brief A specialized Simulation for mio::dabm::Model.
 * @tparam Implementation A class implementing all functions and types marked with the using keyword in Model, see dabm::Model.
 */
template <class Implementation>
class Simulation
{
    using Status = typename Implementation::Status;

public:
    using Model = dabm::Model<Implementation>;

    /**
     * @brief Set up the simulation for a diffusive ABM.
     * @param[in] model An instance of mio::dabm::Model.
     * @param[in] t0 Start time.
     * @param[in] dt Step size of integration.
     */
    Simulation(const Model& model, ScalarType t0 = 0, ScalarType dt = 0.1)
        : m_t0(t0)
        , m_dt(dt)
        , m_model(std::make_unique<Model>(model))
        , m_result(t0, m_model->time_point())
    {
        assert(dt > 0);
        m_current_rates.reserve(m_model->populations.size());
        m_current_events.reserve(m_model->populations.size());
    }

    /**
     * @brief Advance simulation to tmax.
     * This function performs a temporal Gillespie.
     * @param tmax Next stopping point of simulation.
     */
    void advance(const ScalarType t_max)
    {
        // draw time until an adoption takes place
        ScalarType waiting_time = mio::ExponentialDistribution<ScalarType>::get_instance()(m_model->get_rng(), 1.0);
        while (m_t0 < t_max) {
            ScalarType dt             = std::min({m_dt, t_max - m_t0});
            ScalarType remaining_time = dt;
            while (remaining_time > 0) {
                compute_current_rates_and_events(); // lambda_k (aka f-hat(N))
                ScalarType cumulative_adoption_rate =
                    std::accumulate(m_current_rates.begin(), m_current_rates.end(), 0.0); // Lambda
                // status update
                if (waiting_time < cumulative_adoption_rate * remaining_time) {
                    // draw which adoption takes place
                    const size_t event_id =
                        mio::DiscreteDistribution<size_t>::get_instance()(m_model->get_rng(), m_current_rates);
                    Event& event = m_current_events[event_id];
                    // perform adoption
                    m_model->adopt(event.agent, event.new_status);
                    // draw new waiting time
                    remaining_time -= waiting_time / cumulative_adoption_rate;
                    waiting_time = mio::ExponentialDistribution<ScalarType>::get_instance()(m_model->get_rng(), 1.0);
                }
                else {
                    // no event, decrease waiting time
                    waiting_time -= cumulative_adoption_rate * remaining_time;
                    break;
                }
            }
            // position update
            for (auto& agent : m_model->populations) {
                m_model->move(m_t0, dt, agent);
            }
            m_t0 += dt;
            // store result
            m_result.add_time_point(m_t0, m_model->time_point());
        }
    }

    /**
     * @brief Returns the final simulation result.
     * @return A TimeSeries to represent the final simulation result.
     */
    TimeSeries<ScalarType>& get_result()
    {
        return m_result;
    }
    const TimeSeries<ScalarType>& get_result() const
    {
        return m_result;
    }

    /**
     * @brief Returns the model used in the simulation.
     */
    const Model& get_model() const
    {
        return *m_model;
    }
    Model& get_model()
    {
        return *m_model;
    }

private:
    /**
     * @brief Struct defining an adoption event for an agent and target infection state.
     */
    struct Event {
        typename Model::Agent& agent;
        Status new_status;
    };

    /**
     * @brief Calculate current values for m_current_rates and m_current_events.
     */
    inline void compute_current_rates_and_events()
    {
        m_current_rates.clear();
        m_current_events.clear();
        // compute rate for each (agent, status) combination
        for (auto& agent : m_model->populations) {
            for (int s = 0; s < static_cast<int>(Status::Count); s++) {
                Status new_status = static_cast<Status>(s);
                // check if an adoption from the agents status is possible
                auto adoption_rate = m_model->adoption_rate(agent, new_status);
                if (adoption_rate > 0) {
                    // add rate and corresponding event
                    m_current_rates.push_back(adoption_rate);
                    m_current_events.push_back({agent, new_status});
                }
            }
        }
    }

    ScalarType m_t0, m_dt; ///< Start time of the simulation and integration step size.
    std::unique_ptr<Model> m_model; ///< Pointer to the model used in the simulation.
    std::vector<ScalarType> m_current_rates; ///< Current adoption rates.
    std::vector<Event> m_current_events; ///< Contains an event corresponding to each rate in m_current_rates.
    mio::TimeSeries<ScalarType> m_result; ///< Result time series.
};
} //namespace dabm
} // namespace mio

#endif
