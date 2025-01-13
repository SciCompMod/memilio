/* 
* Copyright (C) 2020-2024 German Aerospace Center (DLR-SC)
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

#ifndef MIO_SMM_SIMULATION_H
#define MIO_SMM_SIMULATION_H

#include "memilio/config.h"
#include "smm/model.h"
#include "smm/parameters.h"
#include "memilio/compartments/simulation.h"

namespace mio
{
template <size_t regions, class Status>
class Simulation<ScalarType, smm::Model<regions, Status>>
{
public:
public:
    using Model = smm::Model<regions, Status>;

    /**
     * @brief setup the simulation with an ODE solver
     * @param[in] model: An instance of a compartmental model
     * @param[in] t0 start time
     * @param[in] dt initial step size of integration
     */
    Simulation(Model const& model, ScalarType t0 = 0., ScalarType dt = 1.)
        : m_dt(dt)
        , m_model(std::make_unique<Model>(model))
        , m_result(t0, m_model->get_initial_values())
        , m_internal_time(adoption_rates().size() + transition_rates().size(), t0)
        , m_tp_next_event(adoption_rates().size() + transition_rates().size(), t0)
        , m_waiting_times(adoption_rates().size() + transition_rates().size(), 0)
        , m_current_rates(adoption_rates().size() + transition_rates().size(), 0)
    {
        assert(dt > 0);
        assert(m_waiting_times.size() > 0);
        assert(std::all_of(adoption_rates().begin(), adoption_rates().end(), [](auto&& r) {
            return static_cast<size_t>(r.region) < regions;
        }));
        assert(std::all_of(transition_rates().begin(), transition_rates().end(), [](auto&& r) {
            return static_cast<size_t>(r.from) < regions && static_cast<size_t>(r.to) < regions;
        }));
        // initialize (internal) next event times by random values
        for (size_t i = 0; i < m_tp_next_event.size(); i++) {
            m_tp_next_event[i] += mio::ExponentialDistribution<ScalarType>::get_instance()(m_model->get_rng(), 1.0);
        }
    }

    /**
     * @brief advance simulation to tmax
     * tmax must be greater than get_result().get_last_time_point()
     * @param tmax next stopping point of simulation
     */
    Eigen::Ref<Eigen::VectorXd> advance(ScalarType tmax)
    {
        update_current_rates_and_waiting_times();
        size_t next_event       = determine_next_event(); // index of the next event
        ScalarType current_time = m_result.get_last_time();
        // set in the past to add a new time point immediately
        ScalarType last_result_time = current_time - m_dt;
        // iterate over time
        while (current_time + m_waiting_times[next_event] < tmax) {
            // update time
            current_time += m_waiting_times[next_event];
            // regularily save current state in m_results
            if (current_time > last_result_time + m_dt) {
                last_result_time = current_time;
                m_result.add_time_point(current_time);
                // copy from the previous last value
                m_result.get_last_value() = m_result[m_result.get_num_time_points() - 2];
            }
            // decide event type by index and perform it
            if (next_event < adoption_rates().size()) {
                // perform adoption event
                const auto& rate = adoption_rates()[next_event];
                m_result.get_last_value()[m_model->populations.get_flat_index({rate.region, rate.from})] -= 1;
                m_model->populations[{rate.region, rate.from}] -= 1;
                m_result.get_last_value()[m_model->populations.get_flat_index({rate.region, rate.to})] += 1;
                m_model->populations[{rate.region, rate.to}] += 1;
            }
            else {
                // perform transition event
                const auto& rate = transition_rates()[next_event - adoption_rates().size()];
                m_result.get_last_value()[m_model->populations.get_flat_index({rate.from, rate.status})] -= 1;
                m_model->populations[{rate.from, rate.status}] -= 1;
                m_result.get_last_value()[m_model->populations.get_flat_index({rate.to, rate.status})] += 1;
                m_model->populations[{rate.to, rate.status}] += 1;
            }
            // update internal times
            for (size_t i = 0; i < m_internal_time.size(); i++) {
                m_internal_time[i] += m_current_rates[i] * m_waiting_times[next_event];
            }
            // draw new "next event" time for the occured event
            m_tp_next_event[next_event] +=
                mio::ExponentialDistribution<ScalarType>::get_instance()(m_model->get_rng(), 1.0);
            // precalculate next event
            update_current_rates_and_waiting_times();
            next_event = determine_next_event();
        }
        // copy last result, if no event occurs between last_result_time and tmax
        if (last_result_time < tmax) {
            const auto& val = m_result.get_last_value();
            m_result.add_time_point(tmax, val);
        }
        return m_result.get_last_value();
    }

    /**
     * @brief get_result returns the final simulation result
     * @return a TimeSeries to represent the final simulation result
     */
    TimeSeries<ScalarType>& get_result()
    {
        return m_result;
    }

    /**
     * @brief get_result returns the final simulation result
     * @return a TimeSeries to represent the final simulation result
     */
    const TimeSeries<ScalarType>& get_result() const
    {
        return m_result;
    }

    /**
     * @brief returns the simulation model used in simulation
     */
    const Model& get_model() const
    {
        return *m_model;
    }

    /**
     * @brief returns the simulation model used in simulation
     */
    Model& get_model()
    {
        return *m_model;
    }

private:
    inline constexpr const typename smm::TransitionRates<Status>::Type& transition_rates()
    {
        return m_model->parameters.template get<smm::TransitionRates<Status>>();
    }
    inline constexpr const typename smm::AdoptionRates<Status>::Type& adoption_rates()
    {
        return m_model->parameters.template get<smm::AdoptionRates<Status>>();
    }

    inline void update_current_rates_and_waiting_times()
    {
        size_t i = 0; // shared index for iterating both rates
        for (const auto& rate : adoption_rates()) {
            m_current_rates[i] = m_model->evaluate(rate, m_result.get_last_value());
            m_waiting_times[i] = (m_current_rates[i] > 0)
                                     ? (m_tp_next_event[i] - m_internal_time[i]) / m_current_rates[i]
                                     : std::numeric_limits<ScalarType>::max();
            i++;
        }
        for (const auto& rate : transition_rates()) {
            m_current_rates[i] = m_model->evaluate(rate, m_result.get_last_value());
            m_waiting_times[i] = (m_current_rates[i] > 0)
                                     ? (m_tp_next_event[i] - m_internal_time[i]) / m_current_rates[i]
                                     : std::numeric_limits<ScalarType>::max();
            i++;
        }
    }
    inline size_t determine_next_event()
    {
        return std::distance(m_waiting_times.begin(), std::min_element(m_waiting_times.begin(), m_waiting_times.end()));
    }

    ScalarType m_dt;
    std::unique_ptr<Model> m_model;
    mio::TimeSeries<ScalarType> m_result;

    std::vector<ScalarType> m_internal_time; // internal times of all poisson processes (aka T_k)
    std::vector<ScalarType> m_tp_next_event; // internal time points of next event i after m_internal[i] (aka P_k)
    std::vector<ScalarType> m_waiting_times; // external times between m_internal_time and m_tp_next_event
    std::vector<ScalarType> m_current_rates; // current values of both types of rates
};

} // namespace mio

#endif
