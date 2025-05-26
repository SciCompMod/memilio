/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: David Kerkmann, Sascha Korf, Khoa Nguyen
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
#ifndef MIO_MATH_INTERPOLATION_H_
#define MIO_MATH_INTERPOLATION_H_
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"

#include <cassert>
#include <vector>
#include <algorithm>

namespace mio
{

/**
 * @brief Linear interpolation between two data values.
 *
 * @param[in] x_eval Location to evaluate interpolation.
 * @param[in] x_1 Left node of interpolation.
 * @param[in] x_2 Right node of interpolation.
 * @param[in] y_1 Value at left node.
 * @param[in] y_2 Value at right node.
 * @return Interpolation result.
 */
template <typename X, typename V>
V linear_interpolation(const X& x_eval, const X& x_1, const X& x_2, const V& y1, const V& y2)
{
    const X weight = (x_eval - x_1) / (x_2 - x_1);
    return y1 + weight * (y2 - y1);
}

/**
 * @brief Linear interpolation of a TimeSeries.
 * Assumes that the times in the time series are monotonic increasing. If the times are *strictly* monotonic,
 * this function is continuous in time.
 * If the given interpolation time is outside of the provided time points, this function assumes a constant value of
 * the first/last time point.
 * @param[in] time The time at which to evaluate.
 * @param[in] data Time points to interpolate. At least one is required.
 * @return Interpolation result.
 */
template <class FP>
typename TimeSeries<FP>::Vector linear_interpolation(FP time, const TimeSeries<FP>& data)
{
    assert(data.get_num_time_points() > 0 && "Interpolation requires at least one time point.");
    const auto tp_range = data.get_times();
    // find next time point in data (strictly) after time
    const auto next_tp = std::upper_bound(tp_range.begin(), tp_range.end(), time);
    // interpolate in between values if possible, otherwise return first/last value
    if (next_tp == tp_range.begin()) { // time is before first data point
        return data.get_value(0);
    }
    else if (next_tp == tp_range.end()) { // time is past or equal to last data point
        return data.get_last_value();
    }
    else { // time is in between data points
        const auto i = next_tp - tp_range.begin(); // index of strict upper bound
        return linear_interpolation(time, data.get_time(i - 1), data.get_time(i), data.get_value(i - 1),
                                    data.get_value(i));
    }
}

/**
 * @brief Linear interpolation between two points of a dataset, which is represented by a vector of pairs of node and value.
 * Return 0 if there is less than two points in the dataset.
 * @param[in] vector Vector of pairs of node and value.
 * @param[in] x_eval Location to evaluate interpolation.
 * @return Interpolation result.
 */
template <typename X, typename Y>
Y linear_interpolation_of_data_set(std::vector<std::pair<X, Y>> vector, const X& x_eval)
{
    // If the vector is empty or has only 1 node, throw an error
    if (vector.empty() || vector.size() == 1) {
        log_error("The vector provided in linear_interpolation_of_data_set() must have more than 1 node.");
        return 0.0;
    }

    // Find the corresponding upper position of the node in the data set
    size_t upper_pos = std::upper_bound(vector.begin(), vector.end(), x_eval,
                                        [](double value, const std::pair<X, Y>& node) {
                                            return value <= node.first;
                                        }) -
                       vector.begin();

    if (upper_pos == 0) {
        return vector[upper_pos].second;
    }

    // If the x_eval is between two identifiable nodes in the dataset.
    if (upper_pos < vector.size()) {
        return linear_interpolation(x_eval, vector[upper_pos - 1].first, vector[upper_pos].first,
                                    vector[upper_pos - 1].second, vector[upper_pos].second);
    }
    return 0;
}

} // namespace mio

#endif // MIO_MATH_INTERPOLATION_H_
