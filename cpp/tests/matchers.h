/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#ifndef EPI_TESTS_MATCHERS_H
#define EPI_TESTS_MATCHERS_H

#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/math/floating_point.h"
#include "memilio/io/io.h"
#include "gmock/gmock.h"

/**
 * @brief overload gtest printer function for eigen matrices.
 * @note see https://stackoverflow.com/questions/25146997/teach-google-test-how-to-print-eigen-matrix
 */
template <class M>
struct MatrixPrintWrap : public M {
    friend void PrintTo(const MatrixPrintWrap& m, std::ostream* os)
    {
        if (m.rows() == 1)
        {            
            //print row vector inline
            (*os) << m;
        }
        else if (m.cols() == 1)
        {
            //print col vector inline transposed
            (*os) << m.transpose() << " T";
        }
        else
        {
            //print matrix on its own
            (*os) << '\n' << m;
        }
    }
};

/**
 * @brief wrap m for gtest printing
 * returns a reference to the original object, no copying or moving, mind the lifetime!
 */
template <class M>
const MatrixPrintWrap<M>& print_wrap(const Eigen::EigenBase<M>& m)
{
    return static_cast<const MatrixPrintWrap<M>&>(m);
}

/**
 * gmock matcher, checks if each element of two eigen matrices are within tolerance.
 * @param other matrix to compare
 * @param rtol relative tolerance
 * @param atol absolute tolerance
 * @return matcher that accepts eigen matrix types
 */
MATCHER_P3(MatrixNear, other, rtol, atol,
           "approx. equal to " + testing::PrintToString(print_wrap(other)) + " (rtol = " + testing::PrintToString(rtol) +
               ", atol = " + testing::PrintToString(atol) + ")")
{
    epi::unused(result_listener);
    return ((arg - other).array().abs() <= (atol + rtol * other.array().abs())).all();
}

/**
 * gmock matcher, checks if each element of two eigen matrices are close.
 * @param other matrix to compare
 * @return matcher that accepts eigen matrix types
 */
MATCHER_P(MatrixNear, other,
           "approx. equal to " + testing::PrintToString(print_wrap(other)) + " (rtol = 1e-15, atol = 1e-15)")
{
    epi::unused(result_listener);
    return ((arg - other).array().abs() <= (1e-15 + 1e-15 * other.array().abs())).all();
}

/**
 * gmock matcher, checks if two floating point values are almost equal.
 * @param other value to compare
 * @param rtol relative tolerance
 * @param atol absolute tolerance
 * @return matcher that accepts floating point values.
 */
MATCHER_P3(FloatingPointEqual, other, atol, rtol,
           "approx. equal to " + testing::PrintToString(other) + " (rtol = " + testing::PrintToString(rtol) +
               ", atol = " + testing::PrintToString(atol) + ")")
{
    epi::unused(result_listener);
    return epi::floating_point_equal(arg, other, atol, rtol);
}

/**
 * @brief overload gtest printer function for IOResult.
 * @note see https://stackoverflow.com/questions/25146997/teach-google-test-how-to-print-eigen-matrix
 */
template <class T>
struct IOResultPrintWrap : public epi::IOResult<T> {
    friend void PrintTo(const IOResultPrintWrap& m, std::ostream* os)
    {
        if (m) {
            *os << "Success";
        }
        else {
            *os << "Error: " << m.error().formatted_message();
        }
    }
};

/**
 * @brief wrap an IOResult for gtest printing
 * returns a reference to the original object, no copying or moving, mind the lifetime!
 */
template<class T>
const IOResultPrintWrap<T>& print_wrap(const epi::IOResult<T>& r)
{
    return static_cast<const IOResultPrintWrap<T>&>(r);
}

/**
 * gmock matcher for IOResult.
 * The matcher succeeds if the IOResult represents success.
 * @return matcher that checks an IOResult
 */
MATCHER(IsSuccess, std::string(negation ? "isn't" : "is") + " successful. ")
{
    epi::unused(result_listener);
    return bool(arg);
}

/**
 * gmock matcher that checks whether the elements of a container are linearly spaced.
 * @param b minimum value
 * @param e maximum value
 * @param num_points number of linearly spaced points in [b, e]
 * @return matcher that accepts a stl container
 */
template <class T>
auto ElementsAreLinspace(T b, T e, size_t num_points)
{
    assert(num_points >= 2);

    std::vector<decltype(FloatingPointEqual(std::declval<T>(), std::declval<T>(), std::declval<T>()))> values;
    auto step_size = (e - b) / (num_points - 1);
    for (size_t i = 0; i < num_points; i++) {
        values.push_back(FloatingPointEqual(b + i * step_size, 1e-15 * step_size, 1e-15));
    }
    return testing::ElementsAreArray(values);
}

#endif //EPI_TESTS_MATCHERS_H