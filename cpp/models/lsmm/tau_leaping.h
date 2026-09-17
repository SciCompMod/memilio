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
#ifndef MIO_LSMM_TAU_LEAPING_H
#define MIO_LSMM_TAU_LEAPING_H

#include "lsmm/reconstruction.h"

namespace mio
{
namespace lsmm
{

struct TauLeapResult {
    State state;
    Matrix history;
    Matrix endpoints;
};

/**
 * Approximate bounded binomial tau leaping, with rates frozen at each step's aggregate state.
 * Each individual makes at most one transition per step; all history updates are simultaneous.
 * Counts stay nonnegative and group totals are preserved, including on graphs with cycles.
 * This is not the exact CTMC in Pseudocode S3; assess time-step convergence before using results.
 */
inline TauLeapResult simulate_tau_leaping(const Model& model, const Matrix& initial, RandomNumberGenerator& rng,
                                        double t0, double tmax, double dt)
{
    if (initial.rows() != model.num_states() || initial.cols() == 0 || !std::isfinite(t0) || !std::isfinite(tmax) ||
        tmax < t0 || !std::isfinite(dt) || dt <= 0.) {
        throw std::invalid_argument("Invalid LSMM tau-leaping input.");
    }
    details::validate_counts(initial);
    State z        = initial.rowwise().sum();
    Matrix history = z.asDiagonal();
    const auto& transitions = model.transitions();
    std::vector<std::vector<size_t>> outgoing(static_cast<size_t>(model.num_states()));
    for (size_t r = 0; r < transitions.size(); ++r) {
        outgoing[transitions[r].source].push_back(r);
    }
    std::vector<double> rates(transitions.size());
    using Binomial = DistributionAdapter<std::binomial_distribution<Count>>;
    for (double t = t0; t < tmax;) {
        const double step = std::min(dt, tmax - t);
        if (t + step == t) {
            throw std::invalid_argument("LSMM tau step is too small to advance the current time.");
        }
        for (size_t r = 0; r < transitions.size(); ++r) {
            rates[r] = model.per_capita_rate(r, z);
        }
        Matrix updated = history;
        for (Eigen::Index source = 0; source < model.num_states(); ++source) {
            double total_rate = 0.;
            size_t last       = 0;
            for (auto r : outgoing[source]) {
                total_rate += rates[r];
                if (rates[r] > 0.) {
                    last = r;
                }
            }
            if (!std::isfinite(total_rate)) {
                throw std::domain_error("LSMM total exit rate must be finite.");
            }
            if (total_rate == 0.) {
                continue;
            }
            const double departure_probability = -std::expm1(-total_rate * step);
            for (Eigen::Index origin = 0; origin < model.num_states(); ++origin) {
                Count remaining = Binomial::get_instance()(rng, history(source, origin), departure_probability);
                updated(source, origin) -= remaining;
                double remaining_rate = total_rate;
                for (auto r : outgoing[source]) {
                    if (remaining == 0) {
                        break;
                    }
                    if (rates[r] == 0.) {
                        continue;
                    }
                    const Count count = r == last ? remaining : Binomial::get_instance()(
                                            rng, remaining, std::clamp(rates[r] / remaining_rate, 0., 1.));
                    updated(transitions[r].target, origin) += count;
                    remaining -= count;
                    remaining_rate -= rates[r];
                }
            }
        }
        history = std::move(updated);
        z       = history.rowwise().sum();
        t += step;
    }
    Matrix endpoints = reconstruct(initial, history, rng);
    return {std::move(z), std::move(history), std::move(endpoints)};
}

} // namespace lsmm
} // namespace mio
#endif
