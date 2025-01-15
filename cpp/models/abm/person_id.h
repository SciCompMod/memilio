/* 
* Copyright (C) 2020-2024 MEmilio
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

#ifndef MIO_ABM_PERSON_ID_H
#define MIO_ABM_PERSON_ID_H

#include "memilio/utils/type_safe.h"
#include <limits>

namespace mio
{
namespace abm
{

/// Index for a Person within a Model.
struct MEMILIO_ENABLE_EBO LocalIndex : public mio::TypeSafe<uint32_t, LocalIndex>,
                                       public OperatorComparison<LocalIndex> {
    /// @brief Create an Index.
    LocalIndex(uint32_t index)
        : mio::TypeSafe<uint32_t, LocalIndex>(index)
    {
    }

    /// @brief Create an invalid Index.
    LocalIndex()
        : mio::TypeSafe<uint32_t, LocalIndex>(std::numeric_limits<uint32_t>::max())
    {
    }

    /// @brief Value for invalid Indices.
    const static LocalIndex invalid_index()
    {
        return LocalIndex();
    }
};

/// Unique ID for a Person within a Model.
struct MEMILIO_ENABLE_EBO GlobalID : public mio::TypeSafe<uint64_t, GlobalID>, public OperatorComparison<GlobalID> {
    /// @brief Create an ID.
    GlobalID(uint64_t id)
        : mio::TypeSafe<uint64_t, GlobalID>(id)
    {
    }

    /// @brief Create an invalid ID.
    GlobalID()
        : mio::TypeSafe<uint64_t, GlobalID>(std::numeric_limits<uint64_t>::max())
    {
    }

    /// @brief Value for invalid IDs.
    const static GlobalID invalid_ID()
    {
        return GlobalID();
    }
};

} // namespace abm
} // namespace mio

#endif // MIO_ABM_PERSON_ID_H
