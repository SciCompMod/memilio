/*
* Copyright (C) 2020-2025 MEmilio
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

#include "memilio/config.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/math/floating_point.h"
#include "memilio/io/io.h"
#include "gmock/gmock.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "json/json.h"

namespace Json
{
void PrintTo(const Value& json, std::ostream* os);
} // namespace Json

std::string json_type_to_string(Json::ValueType t);

MATCHER_P(JsonEqual, expected_json, testing::PrintToString(expected_json))
{
    auto match_rec = [&](auto&& match, const Json::Value& a, const Json::Value& b, std::string name) {
        // first check if the types match
        if (a.type() != b.type()) {
            *result_listener << "type mismatch for " << name << ", expected " << json_type_to_string(a.type())
                             << ", actual " << json_type_to_string(b.type());
            return false;
        }
        // handle object types by recursively matching members
        if (a.isObject()) {
            for (auto& key : a.getMemberNames()) {
                if (!b.isMember(key)) {
                    *result_listener << "missing key \"" << key << "\" in " << name;
                    return false;
                }
                if (!match(match, a[key], b[key], name + "[\"" + key + "\"]")) {
                    return false;
                }
            }
        }
        // handle arrays by recursively matching each item
        else if (a.isArray()) {
            if (a.size() != b.size()) {
                *result_listener << "wrong number of items in " << name << ", expected " << a.size() << ", actual "
                                 << b.size();
                return false;
            }
            for (Json::ArrayIndex i = 0; i < a.size(); ++i) {
                if (!match(match, a[i], b[i], name + "[\"" + std::to_string(i) + "\"]")) {
                    return false;
                }
            }
        }
        // handle value types using Json::Value::operator==
        else if (a != b) {
            *result_listener << "value mismatch in " << name << ", expected " << testing::PrintToString(a)
                             << ", actual " << testing::PrintToString(b);
            return false;
        }
        return true;
    };
    return match_rec(match_rec, expected_json, arg, "Json::Value");
}
#endif //MEMILIO_HAS_JSONCPP

/**
 * @brief overload gtest printer function for eigen matrices.
 * @note see https://stackoverflow.com/questions/25146997/teach-google-test-how-to-print-eigen-matrix
 */
template <class M>
struct MatrixPrintWrap : public M {
    friend void PrintTo(const MatrixPrintWrap& m, std::ostream* os)
    {
        if (m.rows() == 1) {
            //print row vector inline
            (*os) << m;
        }
        else if (m.cols() == 1) {
            //print col vector inline transposed
            (*os) << m.transpose() << " T";
        }
        else {
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
           "approx. equal to " + testing::PrintToString(print_wrap(other)) +
               " (rtol = " + testing::PrintToString(rtol) + ", atol = " + testing::PrintToString(atol) + ")")
{
    if (arg.rows() != other.rows() || arg.cols() != other.cols()) {
        *result_listener << "different dimensions";
        return false;
    }
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
    mio::unused(result_listener);
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
    mio::unused(result_listener);
    return mio::floating_point_equal(arg, other, atol, rtol);
}

/**
 * @brief overload gtest printer function for IOResult.
 * @note see https://stackoverflow.com/questions/25146997/teach-google-test-how-to-print-eigen-matrix
 */
template <class T>
struct IOResultPrintWrap : public mio::IOResult<T> {
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
template <class T>
const IOResultPrintWrap<T>& print_wrap(const mio::IOResult<T>& r)
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
    if (arg) {
        return true;
    }

    *result_listener << arg.error().formatted_message();
    return false;
}

/**
 * gmock matcher for IOResult.
 * The matcher succeeds if the IOResult represents failure with the specified status code.
 * @return matcher that checks an IOResult
 */
MATCHER_P(IsFailure, status_code, std::string(negation ? "isn't" : "is") + " failure. ")
{
    if (arg.error().code() == status_code) {
        return true;
    }

    *result_listener << arg.error().formatted_message();
    return false;
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
