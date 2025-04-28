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

#ifndef MIO_ABM_PERSON_ID_H
#define MIO_ABM_PERSON_ID_H

#include "memilio/utils/type_safe.h"
#include <limits>

namespace mio
{
namespace abm
{

/// Unique ID for a Person within a Model.
struct MEMILIO_ENABLE_EBO PersonId : public mio::TypeSafe<uint64_t, PersonId>, public OperatorComparison<PersonId> {
    /// @brief Create an ID.
    PersonId(uint64_t id)
        : mio::TypeSafe<uint64_t, PersonId>(id)
    {
    }

    /// @brief Create an invalid ID.
    PersonId()
        : mio::TypeSafe<uint64_t, PersonId>(std::numeric_limits<uint64_t>::max())
    {
    }

    /// @brief Value for invalid IDs.
    const static PersonId invalid_ID()
    {
        return PersonId();
    }
};

} // namespace abm
} // namespace mio

#endif // MIO_ABM_PERSON_ID_H
