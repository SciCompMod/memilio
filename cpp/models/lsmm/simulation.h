/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#ifndef MIO_LSMM_SIMULATION_H
#define MIO_LSMM_SIMULATION_H

#include "lsmm/model.h"
#include "lsmm/reconstruction.h"

namespace mio
{
namespace lsmm
{

/**
 * Exact autonomous aggregate next-reaction simulation with online history H (Pseudocode S3).
 * RandomNumberGenerator must outlive the simulation. Rates must remain fixed functions of Z.
 * No event sequence or group-resolved paths are stored.
 */
class Simulation
{
public:
    Simulation(const Model& model, const Matrix& initial, RandomNumberGenerator& rng, double t0 = 0.)
        : m_model(model)
        , m_initial(initial)
        , m_rng(rng)
        , m_time(t0)
        , m_rates(model.transitions().size())
        , m_per_capita(model.transitions().size(), std::numeric_limits<double>::quiet_NaN())
        , m_remaining(model.transitions().size())
    {
        if (initial.rows() != model.num_states() || initial.cols() == 0 || !std::isfinite(t0)) {
            throw std::invalid_argument("Invalid LSMM initial state or time.");
        }
        details::validate_counts(initial);
        m_state   = initial.rowwise().sum();
        m_history = m_state.asDiagonal();
        for (auto& clock : m_remaining) {
            clock = ExponentialDistribution<double>::get_instance()(m_rng, 1.);
        }
    }

    /// Advance to tmax, including events at tmax. Subsequent calls retain the residual reaction clocks.
    const State& advance(double tmax)
    {
        if (!std::isfinite(tmax) || tmax < m_time) {
            throw std::invalid_argument("LSMM end time must be finite and not precede the current time.");
        }
        if (!m_rates_initialized && m_time < tmax) {
            for (size_t r = 0; r < m_rates.size(); ++r) {
                update_rate(r, true);
            }
            m_rates_initialized = true;
        }
        while (m_time < tmax) {
            double waiting = std::numeric_limits<double>::infinity();
            size_t next    = 0;
            for (size_t r = 0; r < m_rates.size(); ++r) {
                const double dt = m_rates[r] > 0. ? m_remaining[r] / m_rates[r]
                                                 : std::numeric_limits<double>::infinity();
                if (dt < waiting) {
                    waiting = dt;
                    next    = r;
                }
            }
            const double elapsed = std::min(waiting, tmax - m_time);
            // Remaining integrated propensity P_r - T_r, avoiding growing absolute internal clocks.
            for (size_t r = 0; r < m_rates.size(); ++r) {
                m_remaining[r] = std::max(0., m_remaining[r] - elapsed * m_rates[r]);
            }
            if (waiting > tmax - m_time) {
                m_time = tmax;
                break;
            }
            m_time += waiting;
            const auto& transition = m_model.transitions()[next];
            Count draw = UniformIntDistribution<Count>::get_instance()(m_rng, 0, m_state[transition.source] - 1);
            Eigen::Index origin = 0;
            while (draw >= m_history(transition.source, origin)) {
                draw -= m_history(transition.source, origin++);
            }
            --m_history(transition.source, origin);
            ++m_history(transition.target, origin);
            --m_state[transition.source];
            ++m_state[transition.target];
            m_remaining[next] = ExponentialDistribution<double>::get_instance()(m_rng, 1.);
            ++m_num_events;
            for (const auto& update : m_model.affected_rates(next)) {
                update_rate(update.channel, update.refresh_per_capita);
            }
        }
        return m_state;
    }

    const State& get_state() const
    {
        return m_state;
    }
    const Matrix& get_history() const
    {
        return m_history;
    }
    double get_time() const
    {
        return m_time;
    }
    size_t get_num_events() const
    {
        return m_num_events;
    }
    Matrix sample_endpoints()
    {
        return reconstruct(m_initial, m_history, m_rng);
    }

private:
    void update_rate(size_t channel, bool refresh_per_capita)
    {
        auto& per_capita = m_per_capita[channel];
        if (refresh_per_capita) {
            // Keep an invalidated cache empty until the source pool becomes populated.
            per_capita = std::numeric_limits<double>::quiet_NaN();
        }
        const auto count = m_state[m_model.transitions()[channel].source];
        if (count == 0) {
            m_rates[channel] = 0.;
            return;
        }
        if (std::isnan(per_capita)) {
            per_capita = m_model.per_capita_rate(channel, m_state);
        }
        m_rates[channel] = per_capita * static_cast<double>(count);
        if (!std::isfinite(m_rates[channel])) {
            throw std::domain_error("LSMM rates must be finite and nonnegative.");
        }
    }

    Model m_model;
    Matrix m_initial, m_history;
    State m_state;
    RandomNumberGenerator& m_rng;
    double m_time;
    std::vector<double> m_rates, m_per_capita, m_remaining;
    bool m_rates_initialized = false;
    size_t m_num_events = 0;
};

} // namespace lsmm
} // namespace mio
#endif
