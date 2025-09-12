/* 
* Copyright (C) 2020-2025 MEmilio
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
#ifndef MIO_EPI_ADOPTIONRATE_H
#define MIO_EPI_ADOPTIONRATE_H

#include "memilio/utils/index.h"
#include "memilio/config.h"
#include "memilio/geography/regions.h"

namespace mio
{

/**
 * @brief Struct defining an influence for a second-order adoption.
 * The population having "status" is multiplied with "factor."
 * @tparam Status An infection state enum.
 * @tparam Groups Additional grouping indices.
 */
template <class Status, class... Groups>
struct Influence {
    Status status;
    ScalarType factor;
    mio::regions::Region region;
    std::tuple<Groups...> group_indices{};
};

/**
 * @brief Struct defining a possible status adoption in a Model based on Poisson Processes.
 * The AdoptionRate is considered to be of second-order if there are any "influences".
 * In the d_abm and smm simulations, "from" is implicitly an influence, scaled by "factor". This is multiplied by
 * the sum over all "influences", which scale their "status" with the respective "factor".
 * @tparam Status An infection state enum.
 * @tparam Groups Additional grouping indices.
 */
template <class Status, class... Groups>
struct AdoptionRate {
    Status from; // i
    Status to; // j
    mio::regions::Region region; // k
    ScalarType factor; // gammahat_{ij}^k
    std::vector<Influence<Status, Groups...>> influences; // influences[tau] = ( Psi_{i,j,tau} , gamma_{i,j,tau} )
    std::tuple<Groups...> group_indices{};

    template <class... T>
    AdoptionRate(Status initializer_from, Status initializer_to, mio::regions::Region initializer_region,
                 ScalarType initializer_factor, std::vector<std::tuple<Status, ScalarType>> second_order_influences,
                 T... Ts)
        : from(initializer_from)
        , to(initializer_to)
        , region(initializer_region)
        , factor(initializer_factor)
        , group_indices{}
    {
        for (const auto& [status, scalar] : second_order_influences) {
            influences.emplace_back(Influence<Status, Groups...>{status, scalar, initializer_region});
        }
        mio::unused(Ts...);
    }
    template <class... T>
    AdoptionRate(Status initializer_from, Status initializer_to, mio::regions::Region initializer_region,
                 ScalarType initializer_factor,
                 std::vector<std::tuple<Status, ScalarType, mio::regions::Region, std::tuple<Groups...>>>
                     second_order_influences,
                 T... Ts)
        : from(initializer_from)
        , to(initializer_to)
        , region(initializer_region)
        , factor(initializer_factor)
        , group_indices{}
    {
        for (const auto& [status, scalar, influence_region, groups] : second_order_influences) {
            influences.emplace_back(Influence<Status, Groups...>{status, scalar, influence_region, groups});
        }
        mio::unused(Ts...);
    }
    AdoptionRate(Status initializer_from, Status initializer_to, mio::regions::Region initializer_region,
                 ScalarType initializer_factor, std::tuple<> second_order_influences)
        : from(initializer_from)
        , to(initializer_to)
        , region(initializer_region)
        , factor(initializer_factor)
        , influences{}
        , group_indices{}
    {
        mio::unused(second_order_influences);
    }
};

} // namespace mio

#endif
