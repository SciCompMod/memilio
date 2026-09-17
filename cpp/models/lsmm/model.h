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
#ifndef MIO_LSMM_MODEL_H
#define MIO_LSMM_MODEL_H

#include <Eigen/Core>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mio
{
namespace lsmm
{

using Count  = std::int64_t;
using State  = Eigen::Matrix<Count, Eigen::Dynamic, 1>;
using Matrix = Eigen::Matrix<Count, Eigen::Dynamic, Eigen::Dynamic>;

/// A single-person transition. The nonnegative per-capita rate depends only on the aggregate state.
struct Transition {
    Eigen::Index source, target;
    std::function<double(const State&)> rate;
    /// Complete list of state entries used by rate. Omit for unknown dependencies; empty means constant.
    /// Dependence of the total rate on the source population is included automatically.
    std::optional<std::vector<Eigen::Index>> dependencies = std::nullopt;
};

struct RateUpdate {
    size_t channel;
    bool refresh_per_capita;
};

/// Local aggregate model on an arbitrary finite transition graph, including branches and cycles.
class Model
{
public:
    Model(Eigen::Index states, std::vector<Transition> transitions)
        : m_states(states)
        , m_transitions(std::move(transitions))
    {
        if (states <= 0) {
            throw std::invalid_argument("LSMM needs at least one state.");
        }
        for (const auto& r : m_transitions) {
            if (r.source < 0 || r.source >= states || r.target < 0 || r.target >= states || r.source == r.target ||
                !r.rate) {
                throw std::invalid_argument("Invalid LSMM transition.");
            }
            if (r.dependencies) {
                for (auto state : *r.dependencies) {
                    if (state < 0 || state >= states) {
                        throw std::invalid_argument("Invalid LSMM rate dependency.");
                    }
                }
            }
        }
        for (const auto& event : m_transitions) {
            m_affected_rates.emplace_back();
            const auto changed = [&](Eigen::Index state) {
                return state == event.source || state == event.target;
            };
            for (size_t channel = 0; channel < m_transitions.size(); ++channel) {
                const auto& r = m_transitions[channel];
                const bool refresh = !r.dependencies ||
                                     std::any_of(r.dependencies->begin(), r.dependencies->end(), changed);
                if (refresh || changed(r.source)) {
                    m_affected_rates.back().push_back({channel, refresh});
                }
            }
        }
    }

    Eigen::Index num_states() const
    {
        return m_states;
    }

    const std::vector<Transition>& transitions() const
    {
        return m_transitions;
    }

    /// Rate updates required after the given channel fires, including changes of source counts.
    const std::vector<RateUpdate>& affected_rates(size_t channel) const
    {
        return m_affected_rates[channel];
    }

    double per_capita_rate(size_t channel, const State& z) const
    {
        const auto& r = m_transitions[channel];
        // Empty source pools need no rate evaluation (e.g. avoid an empty-population division).
        const double value = z[r.source] == 0 ? 0. : r.rate(z);
        if (!std::isfinite(value) || value < 0. || !std::isfinite(value * static_cast<double>(z[r.source]))) {
            throw std::domain_error("LSMM rates must be finite and nonnegative.");
        }
        return value;
    }

private:
    Eigen::Index m_states;
    std::vector<Transition> m_transitions;
    std::vector<std::vector<RateUpdate>> m_affected_rates;
};

namespace details
{
inline void validate_counts(const Matrix& counts)
{
    // This bound also keeps integer products in the hypergeometric mode calculation within int64_t.
    Count total = 0;
    for (Eigen::Index i = 0; i < counts.size(); ++i) {
        const auto value = counts.data()[i];
        if (value < 0 || value > std::numeric_limits<int>::max() - total) {
            throw std::invalid_argument("LSMM counts must be nonnegative with total at most INT_MAX.");
        }
        total += value;
    }
}
} // namespace details

} // namespace lsmm
} // namespace mio
#endif
