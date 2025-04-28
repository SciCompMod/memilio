/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding
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

#ifndef MIO_ABM_LOCATION_ID_H
#define MIO_ABM_LOCATION_ID_H

#include "memilio/utils/type_safe.h"
#include <limits>

namespace mio
{
namespace abm
{

/// Unique identifier for a Location within a Model.
struct MEMILIO_ENABLE_EBO LocationId : public mio::TypeSafe<uint32_t, LocationId>,
                                       public OperatorComparison<LocationId> {
    /// @brief Create an ID.
    LocationId(uint32_t id)
        : mio::TypeSafe<uint32_t, LocationId>(id)
    {
    }

    /// @brief Create an invalid ID.
    LocationId()
        : mio::TypeSafe<uint32_t, LocationId>(std::numeric_limits<uint32_t>::max())
    {
    }

    /// @brief Value for invalid IDs.
    const static LocationId invalid_id()
    {
        return LocationId();
    }
};

} // namespace abm
} // namespace mio

#endif // MIO_ABM_LOCATION_ID_H
