/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Carlotta Gerstein
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
#ifndef ODESEIRMETAPOP_REGIONS_H
#define ODESEIRMETAPOP_REGIONS_H

#include "memilio/utils/index.h"

namespace mio
{
namespace oseirmetapop
{

/**
 * @brief The Region struct is used as a dynamically
 * sized tag for all region dependent categories
 */
struct Region : public Index<Region> {
    Region(size_t val)
        : Index<Region>(val)
    {
    }
};

} // namespace oseirmetapop
} // namespace mio

#endif
