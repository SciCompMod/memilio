/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele
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
#ifndef MIO_MATH_FLOATING_POINT_H
#define MIO_MATH_FLOATING_POINT_H

#include <algorithm>
#include <cmath>
#include <limits>

namespace mio
{

/**
 * maximum absolute value of two numbers.
 * @param v1 first number
 * @param v2 second number
 * @return maximum absolute value between v1 and v2
 */
template <class T>
T abs_max(T v1, T v2)
{
    using std::abs;
    using std::max;
    return max<T>(abs(v1), abs(v2));
}

/**
 * compare two floating point values for equality with tolerances.
 * Use absolute tolerance for comparisons with zero or if you know the magnitude of the values.
 * Otherwise use relative tolerance. If unsure, use both.
 * @param v1 first floating point value
 * @param v2 second floating point value
 * @param abs_tol maximum allowed absolute difference, default 0.
 * @param rel_tol maximum allowed relative difference, default numeric_limits::min.
 * @return true if v1 is within the specified relative OR absolute tolerance of v2  
 */
template <class T>
bool floating_point_equal(T v1, T v2, T abs_tol = 0, T rel_tol = std::numeric_limits<T>::min())
{
    using std::abs;
    auto diff = abs(v1 - v2);
    return diff <= abs_tol || diff <= abs_max(v1, v2) * rel_tol;
}

/**
 * compare two floating point values with tolerances.
 * v1 < v2 if 
 *  a) v1 not == v2 within tolerances and 
 *  b) v1 not > v2.
 * Use absolute tolerance for comparisons with zero or if you know the magnitude of the values.
 * Use relative tolerance (or both) otherwise.
 * @param v1 first floating point value
 * @param v2 second floating point value
 * @param abs_tol maximum allowed absolute difference for equality, default 0. 
 * @param rel_tol maximum allowed relative difference for equality, default numeric_limits::min.
 * @return true if v1 is less than v2 and not within relative or absolute tolerance of v2.
 */
template <class T>
bool floating_point_less(T v1, T v2, T abs_tol = 0, T rel_tol = std::numeric_limits<T>::min())
{
    auto diff = v1 - v2;
    return diff < -(abs_tol + rel_tol * abs_max(v1, v2));
}

/**
 * compare two floating point values with tolerances.
 * v1 > v2 if 
 *  a) v1 not == v2 within tolerances AND 
 *  b) v1 not < v2.
 * Use absolute tolerance for comparisons with zero or if you know the magnitude of the values.
 * Use relative tolerance (or both) otherwise.
 * @param v1 first floating point value
 * @param v2 second floating point value
 * @param abs_tol maximum allowed absolute difference, default 0. 
 * @param rel_tol maximum allowed relative difference, default numeric_limits::min.
 * @return true if v1 is greater than v2 and not within absolute or relative tolerance of v2.
 */
template <class T>
bool floating_point_greater(T v1, T v2, T abs_tol = 0, T rel_tol = std::numeric_limits<T>::min())
{
    return floating_point_less(v2, v1, abs_tol, rel_tol);
}

/**
 * compare two floating point values with tolerances.
 * v1 <= v2 if 
 *  a) v1 < v2 OR 
 *  b) v1 == v2 within tolerances.
 * Use absolute tolerance for comparisons with zero or if you know the magnitude of the values.
 * Use relative tolerance (or both) otherwise.
 * @param v1 first floating point value
 * @param v2 second floating point value
 * @param abs_tol maximum allowed absolute difference, default 0. 
 * @param rel_tol maximum allowed relative difference, default numeric_limits::min.
 * @return true if v1 is less than v2 or within relative or absolute tolerances of v2.
 */
template <class T>
bool floating_point_less_equal(T v1, T v2, T abs_tol = 0, T rel_tol = std::numeric_limits<T>::min())
{
    return !floating_point_greater(v1, v2, abs_tol, rel_tol);
}

/**
 * compare two floating point values with tolerances.
 * v1 >= v2 if 
 *  a) v1 > v2 OR 
 *  b) v1 == v2 within tolerances.
 * Use absolute tolerance for comparisons with zero or if you know the magnitude of the values.
 * Use relative tolerance (or both) otherwise.
 * @param v1 first floating point value
 * @param v2 second floating point value
 * @param abs_tol maximum allowed absolute difference, default 0. 
 * @param rel_tol maximum allowed relative difference, default numeric_limits::min.
 * @return true if v1 is greater than v2 or within absolute or relative tolerance of v2.
 */
template <class T>
bool floating_point_greater_equal(T v1, T v2, T abs_tol = 0, T rel_tol = std::numeric_limits<T>::min())
{
    return !floating_point_less(v1, v2, abs_tol, rel_tol);
}

/**
 * @brief Rounds a value to the nearest nth decimal place.
 * 
 * @param x The double value to be rounded.
 * @param n The number of decimal places we want to round to.
 * @return The rounded double value to n decimal digits.
 */
template <class T>
T round_nth_decimal(T x, size_t n)
{
    using std::round;
    T factor = std::pow(10.0, n);
    return round(x * factor) / factor;
}

} // namespace mio

#endif // MIO_MATH_FLOATING_POINT_H
