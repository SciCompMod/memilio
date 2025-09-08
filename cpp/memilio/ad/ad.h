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
