/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: René Schmieding, Julia Bicker
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

#ifndef MIO_SMM_PARAMETERS_H
#define MIO_SMM_PARAMETERS_H

#include "memilio/config.h"
#include "memilio/geography/regions.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/epidemiology/adoption_rate.h"

namespace mio
{
namespace smm
{

/**
 * @brief A vector of AdoptionRate%s, see mio::AdoptionRate
 * @tparam FP A floating point type, e.g., double.
 * @tparam Status A MultiIndex, containing the infection state enum and all stratification groups.
 * @tparam Region A MultiIndex for spatial stratification.
 */
template <typename FP, class Status, class Region>
struct AdoptionRates {
    using Type = std::vector<AdoptionRate<FP, Status, Region>>;
    const static std::string name()
    {
        return "AdoptionRates";
    }
};

/**
 * @brief Struct defining a possible regional transition in a Model based on Poisson Processes.
 * @tparam FP A floating point type, e.g., double.
 * @tparam Status A MultiIndex, containing the infection state enum and all stratification groups.
 * @tparam Region A MultiIndex for spatial stratification.
 */
template <typename FP, class Status, class Region = mio::regions::Region>
struct TransitionRate {
    Status status; // i
    Region from; // k
    Region to; // l
    FP factor; // lambda_i^{kl}
};

/**
 * @brief A vector of TransitionRate%s, see mio::TransitionRate
 * @tparam FP A floating point type, e.g., double.
 * @tparam Status A MultiIndex, containing the infection state enum.
 * @tparam Region A MultiIndex for spatial stratification.
 */
template <typename FP, class Status, class Region>
struct TransitionRates {
    using Type = std::vector<TransitionRate<FP, Status, Region>>;
    const static std::string name()
    {
        return "TransitionRates";
    }
};

template <typename FP, class Status, class Region>
using ParametersBase = mio::ParameterSet<AdoptionRates<FP, Status, Region>, TransitionRates<FP, Status, Region>>;

} // namespace smm

} // namespace mio

#endif
