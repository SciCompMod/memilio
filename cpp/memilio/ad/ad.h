/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Julian Litz
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

#ifndef MIO_AD_H
#define MIO_AD_H

#include "ad/ad.hpp"
#include "memilio/math/eigen.h"

#include <cmath>
#include <limits>

// Extend automatic differentiation (AD) library to support std::round.
namespace ad
{
namespace internal
{
using std::round;
template <class AD_TAPE_REAL, class DATA_HANDLER_1>
static inline double round(const ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER_1>& x)
{
    return round(x._value());
}
template <class AD_TAPE_REAL, class A1_T1, class A1_T2, class A1_OP>
static inline double round(const ad::internal::binary_intermediate_aa<AD_TAPE_REAL, A1_T1, A1_T2, A1_OP>& x)
{
    return round(x._value());
}
template <class AD_TAPE_REAL, class A1_T1, class A1_OP>
static inline double round(const ad::internal::binary_intermediate_ap<AD_TAPE_REAL, A1_T1, A1_OP>& x)
{
    return round(x._value());
}
template <class AD_TAPE_REAL, class A1_T2, class A1_OP>
static inline double round(const ad::internal::binary_intermediate_pa<AD_TAPE_REAL, A1_T2, A1_OP>& x)
{
    return round(x._value());
}
template <class AD_TAPE_REAL, class A1_T, class A1_OP>
static inline double round(const ad::internal::unary_intermediate<AD_TAPE_REAL, A1_T, A1_OP>& x)
{
    return round(x._value());
}
} // namespace internal
} // namespace ad

// Allow std::numeric_limits to work with AD types.
template <class FP, class DataHandler>
struct std::numeric_limits<ad::internal::active_type<FP, DataHandler>> : public numeric_limits<FP> {
};

// Ensures that Eigen recognizes your AD types as valid scalars.
namespace Eigen
{
template <class FP, class DataHandler>
struct NumTraits<ad::internal::active_type<FP, DataHandler>>
    : GenericNumTraits<ad::internal::active_type<FP, DataHandler>> {
    using Scalar     = ad::internal::active_type<FP, DataHandler>;
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

#endif // MIO_AD_H
