/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Rene Schmieding
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
#ifndef MIO_MATH_MATH_UTILS_H
#define MIO_MATH_MATH_UTILS_H

#include "memilio/config.h"
#include "memilio/io/io.h"
#include "memilio/math/eigen.h" // IWYU pragma: keep

namespace mio
{

/**
 * @brief Map a vector onto nonnegative values while preserving its nonnegative sum.
 * Requires the vector to have a nonnegative sum to function correctly.
 * In case of a negative sum, this function distributes the absolute sum evenly over all entries and returns an error.
 * @param[in,out] x The vector to map.
 * @param[in] tolerance Absolute tolerance to check if a value is 0.
 * @return Returns an error if the sum was negative. 
 */
template <typename FP>
IOResult<void> map_to_nonnegative(Eigen::Ref<Eigen::VectorX<FP>> x, const FP tolerance = Limits<FP>::zero_tolerance())
{
    assert(tolerance > 0);
    FP positive{0.0}, negative{0.0};
    for (auto& v : x) {
        if (v >= FP{0.0}) {
            positive += v;
        }
        else {
            negative -= v;
            v = FP{0.0};
        }
    }
    // check if there is anything to do
    if (negative > tolerance) {
        // check if sum >= 0. otherwise, we cannot preserve it
        if (positive >= negative) [[likely]] {
            x = x * (1 - negative / positive);
        }
        else {
            x.array() = (negative - positive) / x.size(); // = abs(positive - negative) / x.size()
            return failure(StatusCode::InvalidValue,
                           "Failed to map vector to nonnegative values, total sum is negative.");
        }
    }
    return success();
}

/**
 * @brief Evaluate an intermediate expression to its underlying type, if necessary.
 * @param x An intermediate expression, resulting e.g. from `auto x = y + z`.
 * @tparam Type Underlying type of the expression, e.g. an AD type. Usually this is set to `FP`.
 * @tparam Intermediate Type of the expression result. Do not specify this template, let the compiler deduce it.
 * @return Either casts the intermediate expression to Type, or forwards it.
 *
 * The main purpose of this function is reconcile the handling of intermediate types from Eigen and AD.
 * Eigen wants to return intermediates for optimization, while AD values must be evaluated (due to using references).
 */
template <typename Type, typename Intermediate>
inline auto evaluate_intermediate(Intermediate&& x)
{
    if constexpr (is_ad_type_v<Type>) {
        return static_cast<Type>(x);
    }
    else {
        return std::forward<Intermediate>(x);
    }
}

} // namespace mio

#endif // MIO_MATH_MATH_UTILS_H
