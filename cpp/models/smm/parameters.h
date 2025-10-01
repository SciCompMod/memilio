/*
* Copyright (C) 2020-2025 German Aerospace Center (DLR-SC)
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
#include <tuple>

namespace mio
{
namespace smm
{

/**
 * @brief A vector of AdoptionRate%s, see mio::AdoptionRate
 * @tparam Status An infection state enum.
 */
template <typename FP, class Status, class... Groups>
struct AdoptionRates {
    using Type = std::vector<AdoptionRate<FP, Status, Groups...>>;
    const static std::string name()
    {
        return "AdoptionRates";
    }
};

/**
 * @brief Struct defining a possible regional transition in a Model based on Poisson Processes.
 * @tparam Status An infection state enum.
 */
template <typename FP, class Status, class... Groups>
struct TransitionRate {
    Status status; // i
    mio::regions::Region from; // k
    mio::regions::Region to; // l
    FP factor; // lambda_i^{kl}
    std::tuple<Groups...> group_indices_from{};
    std::tuple<Groups...> group_indices_to{};
};
template <typename FP, class Status, class... Groups>
struct TransitionRates {
    using Type = std::vector<TransitionRate<FP, Status, Groups...>>;
    const static std::string name()
    {
        return "TransitionRates";
    }
};

template <typename FP, class Status, class... Groups>
using ParametersBase = mio::ParameterSet<AdoptionRates<FP, Status, Groups...>, TransitionRates<FP, Status, Groups...>>;

} // namespace smm

} // namespace mio

#endif
