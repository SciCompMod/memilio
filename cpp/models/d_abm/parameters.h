/* 
* Copyright (C) 2020-2024 German Aerospace Center (DLR-SC)
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

#ifndef MIO_d_ABM_PARAMETERS_H
#define MIO_d_ABM_PARAMETERS_H

#include "memilio/config.h"
#include "memilio/utils/index.h"

namespace mio
{
namespace dabm
{
/// @biref Index for enumerating subsections (cities, counties, etc.) of the modelled area.
struct Region : public mio::Index<Region> {
    Region(const size_t num_regions)
        : mio::Index<Region>(num_regions)
    {
    }
};

/**
 * @brief Struct defining a possible status adoption in a Model based on Poisson Processes.
 * The AdoptionRate is considered to be of second order, if there are any "influences" with corresponding "factors".
 * In the d_abm and smm simulations, "from" is implicitly an influence, scaled by "factor". This is multiplied by
 * the sum over all "influences" scaled by their respective "factors".
 * @tparam Status An infection state enum.
 */
template <class Status>
struct AdoptionRate {
    Status from; // i
    Status to; // j
    Region region; // k
    ScalarType factor; // gammahat_{ij}^k
    std::vector<Status> influences;
    std::vector<ScalarType> factors;
};

template <class Status>
struct AdoptionRates {
    using Type = std::vector<AdoptionRate<Status>>;
    const static std::string name()
    {
        return "AdoptionRates";
    }
};

} // namespace dabm

} // namespace mio

#endif
