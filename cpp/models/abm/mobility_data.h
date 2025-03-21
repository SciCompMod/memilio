/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Sascha Korf
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

#ifndef ABM_MOBILITY_DATA_H
#define ABM_MOBILITY_DATA_H

#include <cstdint>

namespace mio
{
namespace abm
{

/**
 * @brief Mode of Transport.
 */
enum class TransportMode : uint32_t
{
    Bike,
    CarDriver,
    CarPassenger,
    PublicTransport,
    Walking,
    Other,
    Unknown,
    Count //last!!
};

/**
 * @brief Type of the activity.
 */
enum class ActivityType : uint32_t
{
    Workplace,
    Education,
    Shopping,
    Leisure,
    PrivateMatters,
    OtherActivity,
    Home,
    UnknownActivity
};

} // namespace abm
} // namespace mio

#endif //ABM_MOBILITY_DATA_H
