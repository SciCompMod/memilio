/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Julia Bicker
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
#ifndef MIO_ABM_ACTIVITY_TYPE_H
#define MIO_ABM_ACTIVITY_TYPE_H

#include <cstdint>

namespace mio
{
namespace abm
{

/**
 * @brief Type of an Activity. This is used to determine the type of an Activity that a Person does at a Location. It is similar to LocationType, but is not necessarily the same as Persons can do different activities at the same location e.g. "Work" and "School" at a Location of LocationType "School".
 */
enum class ActivityType : std::uint32_t
{
    Home = 0,
    School,
    Work,
    Recreation, // TODO: differentiate different kinds
    BasicsShop, // groceries and other necessities
    Hospital,
    ICU,
    Cemetery, // Location for all the dead persons. It is created once for the Model.
    PublicTransport,

    Count, //last!
    Invalid
};

} // namespace abm
} // namespace mio

#endif
