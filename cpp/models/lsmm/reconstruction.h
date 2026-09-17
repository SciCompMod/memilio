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
#ifndef MIO_LSMM_RECONSTRUCTION_H
#define MIO_LSMM_RECONSTRUCTION_H

#include "lsmm/model.h"
#include "memilio/utils/random_number_generator.h"
#include <algorithm>

namespace mio
{
namespace lsmm
{
namespace details
{
inline Count hypergeometric(Count population, Count successes, Count draws, RandomNumberGenerator& rng)
{
    if (draws == 0 || successes == 0) {
        return 0;
    }
    if (draws == population || successes == population) {
        return std::min(draws, successes);
    }
    const Count lower = std::max(Count(0), draws - (population - successes));
    const Count upper = std::min(draws, successes);
    const Count mode  = std::clamp((draws + 1) * (successes + 1) / (population + 2), lower, upper);
    // Adjacent probability ratios; relative mass at the mode is one (no factorials).
    const auto right = [=](Count k) -> long double {
        return static_cast<long double>(successes - k) * (draws - k) /
               (static_cast<long double>(k + 1) * (population - successes - draws + k + 1));
    };
    const auto left = [=](Count k) -> long double {
        return static_cast<long double>(k) * (population - successes - draws + k) /
               (static_cast<long double>(successes - k + 1) * (draws - k + 1));
    };
    // Moving away from the mode, ratios decrease. The geometric bound limits omitted tail mass
    // to long-double precision. This avoids scanning a population-sized zero-probability tail.
    long double mass = 1., weight = 1.;
    Count hi = mode, lo = mode;
    const auto epsilon = std::numeric_limits<long double>::epsilon();
    while (hi < upper) {
        weight *= right(hi++);
        mass += weight;
        if (weight <= epsilon * mass * (1 - right(hi))) {
            break;
        }
    }
    weight = 1.;
    while (lo > lower) {
        weight *= left(lo--);
        mass += weight;
        if (weight <= epsilon * mass * (1 - left(lo))) {
            break;
        }
    }
    long double draw = UniformDistribution<double>::get_instance()(rng, 0., 1.) * mass;
    if (draw < 1.) {
        return mode;
    }
    draw -= 1.;
    weight = 1.;
    for (Count k = mode; k < hi;) {
        weight *= right(k++);
        if (draw < weight) {
            return k;
        }
        draw -= weight;
    }
    weight = 1.;
    for (Count k = mode; k > lo;) {
        weight *= left(k--);
        if (draw < weight) {
            return k;
        }
        draw -= weight;
    }
    return lo; // Last retained category receives the floating-point residual.
}
} // namespace details

/**
 * Direct endpoint sampler from Pseudocode S3.
 * initial(state, group) stores initial counts; history(current state, initial state) stores H.
 * The resulting columns conserve group totals and sum to the row sums of H.
 */
inline Matrix reconstruct(const Matrix& initial, const Matrix& history, RandomNumberGenerator& rng)
{
    if (initial.rows() == 0 || initial.cols() == 0 || history.rows() != initial.rows() ||
        history.cols() != initial.rows()) {
        throw std::invalid_argument("Invalid LSMM reconstruction dimensions.");
    }
    details::validate_counts(initial);
    details::validate_counts(history);
    if (history.colwise().sum().transpose() != initial.rowwise().sum()) {
        throw std::invalid_argument("History column totals must equal initial aggregate counts.");
    }
    Matrix endpoints = Matrix::Zero(initial.rows(), initial.cols());
    for (Eigen::Index j = 0; j < initial.rows(); ++j) {
        State residual = history.col(j);
        for (Eigen::Index g = 0; g + 1 < initial.cols(); ++g) {
            Count population = residual.sum();
            Count draws      = initial(j, g);
            for (Eigen::Index i = 0; i < initial.rows() && draws > 0; ++i) {
                const Count count = details::hypergeometric(population, residual[i], draws, rng);
                population -= residual[i];
                residual[i] -= count;
                draws -= count;
                endpoints(i, g) += count;
            }
        }
        endpoints.col(initial.cols() - 1) += residual;
    }
    return endpoints;
}

} // namespace lsmm
} // namespace mio
#endif
