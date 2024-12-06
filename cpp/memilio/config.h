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

/**
 * Configuration of memilio library.
 */

#ifndef MIO_CONFIG_H
#define MIO_CONFIG_H

#include "memilio/config_internal.h"
#include <type_traits>
#include "ad/ad.hpp"

using ScalarType = double;

namespace mio
{
/**
 * @brief Type specific limits for floating-points.
 *
 * Specializations are provided for `float` and `double` types. In order to use
 * other floating point types, a new specialization for that type must be added.
 * 
 * Member functions of the unspecialized class are disabled. Using them results
 * in compile time errors.
 * @{
 */
/**
 * @tparam FP A floating point type.
 */
template <typename FP>
struct Limits {
    static constexpr FP zero_tolerance() = delete;
};

template <>
struct Limits<float> {
    /// @brief Returns the limit under which a float may be rounded down to zero.
    static constexpr float zero_tolerance()
    {
        return 1e-6f;
    }
};

template <>
struct Limits<double> {
    /// @brief Returns the limit under which a double may be rounded down to zero.
    static constexpr double zero_tolerance()
    {
        return 1e-12;
    }
};

template <class AD_TAPE_REAL, class DATA_HANDLER>
struct Limits<ad::internal::active_type<AD_TAPE_REAL, DATA_HANDLER>> {
    /// @brief Returns the limit under which an ad::internal::active_type may be rounded down to zero.
    static constexpr AD_TAPE_REAL zero_tolerance()
    {
        return Limits<AD_TAPE_REAL>::zero_tolerance();
    }
};
/** @} */

} // namespace mio
#endif
