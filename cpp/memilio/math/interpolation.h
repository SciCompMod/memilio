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
#include <iostream>
#include <vector>

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
 * @param[in] node The node whose value need to be found.
 * @param[out] unnamed Interpolation result.
 */
template <typename X, typename Y>
Y linear_interpolation_of_data_set(const std::vector<std::pair<X, Y>>& vector, const X& node)
{
    // If the vector is empty or has only 1 node, return 0
    if (vector.empty() || vector.size() == 1) {
        return 0.0;
    }

    std::vector<std::pair<X, X>> copy_vector(vector);
    sort(copy_vector.begin(), copy_vector.end());

    // Find the corresponding section of the node in the data set
    size_t counter = 0;
    while ((counter < copy_vector.size() - 1) && (copy_vector[counter].first < node)) {
        counter++;
    }
    std::cout << "Counter: " << counter << "\n";
    std::cout << "Node: " << node << "\n";
    // If the current days interval are between two identifiable points in the dataset.
    if (node <= copy_vector[counter].first && counter > 0) {
        return copy_vector[counter - 1].second + (copy_vector[counter - 1].second - copy_vector[counter].second) /
                                                     (copy_vector[counter - 1].first - copy_vector[counter].first) *
                                                     (node - copy_vector[counter - 1].first);
    }
    return 0.0;
}

} // namespace mio

#endif
