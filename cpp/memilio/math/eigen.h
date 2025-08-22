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
#ifndef MIO_UTILS_EIGEN_H
#define MIO_UTILS_EIGEN_H

// this file wraps includes from eigen3 library to disable warnings
// Eigen is treated as a system library through cmake since #1168, so ALL warnings within the library are disabled

#include <Eigen/Core>

#include "ad/ad.hpp"

// The following Eigen::NumTraits specializations are required so that Eigen
// recognizes ad::gt1s<double>::type and ad::ga1s<double>::type as valid scalar types.
// Without these, Eigen's internal mechanisms (such as .abs(), .cwiseAbs(), etc.)
// will not work correctly with AD types, because Eigen queries NumTraits to determine
// type properties and which mathematical functions are available.
// This enables Eigen to use your AD types in arrays, matrices, and all vectorized math.

namespace Eigen
{
template <>
struct NumTraits<ad::gt1s<double>::type> : GenericNumTraits<ad::gt1s<double>::type> {
    using Scalar     = ad::gt1s<double>::type;
    using Real       = Scalar;
    using NonInteger = Scalar;
    using Nested     = Scalar;
    enum
    {
        IsComplex             = 0,
        IsInteger             = 0,
        IsSigned              = 1,
        RequireInitialization = 1,
        ReadCost              = 1,
        AddCost               = 3,
        MulCost               = 3
    };
};

template <>
struct NumTraits<ad::ga1s<double>::type> : GenericNumTraits<ad::ga1s<double>::type> {
    using Scalar     = ad::ga1s<double>::type;
    using Real       = Scalar;
    using NonInteger = Scalar;
    using Nested     = Scalar;
    enum
    {
        IsComplex             = 0,
        IsInteger             = 0,
        IsSigned              = 1,
        RequireInitialization = 1,
        ReadCost              = 1,
        AddCost               = 3,
        MulCost               = 3
    };
};

} // namespace Eigen

#endif // MIO_UTILS_EIGEN_H
