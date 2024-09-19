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

using ScalarType = double;

template <typename FP>
struct ZeroTolFP {
    static_assert(std::is_floating_point<FP>::value, "ZeroTolFP is only valid for floating point types.");
    static constexpr FP value = FP(1e-6);
};

template <>
struct ZeroTolFP<float> {
    static constexpr float value = 1e-6f;
};

template <>
struct ZeroTolFP<double> {
    static constexpr double value = 1e-12;
};

#endif
