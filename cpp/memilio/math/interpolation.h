/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: David Kerkmann, Sascha Korf
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

} // namespace mio

#endif
