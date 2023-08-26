/*
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_
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
 * @param[out] unnamed Interpolation result.
 */
template <typename X, typename V>
auto linear_interpolation(const X& x_eval, const X& x_1, const X& x_2, const V& y1, const V& y2)
{
    auto weight = (x_eval - x_1) / (x_2 - x_1);
    return y1 + weight * (y2 - y1);
}

/**
 * @brief Linear interpolation between two point of a dataset, which is represented by a vector of pairs of node and value.
 * @param[in] vector Vector of pairs of node and value.
 * @param[in] x_eval Location to evaluate interpolation.
 * @param[out] unnamed Interpolation result.
 */
template <typename X, typename Y>
Y linear_interpolation_of_data_set(const std::vector<std::pair<X, Y>> vector, const X& x_eval)
{
    // If the vector is empty or has only 1 node, return 0
    if (vector.empty() || vector.size() == 1) {
        throw std::invalid_argument(
            "The vector provided in linear_interpolation_of_data_set() must have larger than 2 nodes");
        return 0.0;
    }

    std::vector<std::pair<X, Y>> copy_vector(vector);
    sort(copy_vector.begin(), copy_vector.end());

    // Find the corresponding upper position of the node in the data set
    size_t upper_pos = std::upper_bound(copy_vector.begin(), copy_vector.end(), x_eval,
                                        [](double value, const std::pair<X, Y>& node) {
                                            return value <= node.first;
                                        }) -
                       copy_vector.begin();

    // If the x_eval are between two identifiable nodes in the dataset.
    if (upper_pos < copy_vector.size() && upper_pos > 0) {
        return linear_interpolation(x_eval, copy_vector[upper_pos - 1].first, copy_vector[upper_pos].first,
                                    copy_vector[upper_pos - 1].second, copy_vector[upper_pos].second);
    }
    return 0.0;
}

} // namespace mio

#endif
