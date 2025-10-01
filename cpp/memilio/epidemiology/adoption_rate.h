/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: René Schmieding, Julia Bicker, Kilian Volmer
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
#include <limits>
#include <tuple>
#include <optional>

namespace mio
{

/**
 * @brief Struct defining a possible status adoption in a Model based on Poisson Processes.
 * The AdoptionRate is considered to be of second-order if there are any "influences".
 * In the d_abm and smm simulations, "from" is implicitly an influence, scaled by "factor". This is multiplied by
 * the sum over all "influences", which scale their "status" with the respective "factor".
 * @tparam Status An infection state enum.
 * @tparam Groups Additional grouping indices.
 */
template <typename FP, class Status, class... Groups>
struct AdoptionRate {

    /**
 * @brief Struct defining an influence for a second-order adoption.
 * The population having "status" is multiplied with "factor."
 * @tparam status An infection state enum.
 * @tparam factor Scaling factor for the influence.
 * @tparam 
 * @tparam Groups Additional grouping indices.
 */
    struct Influence {
        Status status;
        FP factor;
        std::optional<mio::regions::Region> region = std::nullopt;
        std::tuple<Groups...> group_indices{};
    };

    Status from; // i
    Status to; // j
    mio::regions::Region region; // k
    FP factor; // gammahat_{ij}^k
    std::vector<Influence> influences; // influences[tau] = ( Psi_{i,j,tau} , gamma_{i,j,tau} )
    std::tuple<Groups...> group_indices{};
};

} // namespace mio

#endif
