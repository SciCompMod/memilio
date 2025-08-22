/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, Khoa Nguyen, Sascha Korf, Carlotta Gerstein
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
#ifndef MIO_ABM_LOCATION_TYPE_H
#define MIO_ABM_LOCATION_TYPE_H

#include <cstdint>

namespace mio
{
namespace abm
{

/**
 * @brief Type of a Location.
 */
enum class LocationType : std::uint32_t
{
    Home = 0,
    School,
    Work,
    SocialEvent, // TODO: differentiate different kinds
    BasicsShop, // groceries and other necessities
    Hospital,
    ICU,
    Car,
    PublicTransport,
    TransportWithoutContact, // all ways of travel with no contact to other people, e.g. biking or walking
    Cemetery, // Location for all the dead persons. It is created once for the Model.

    Count, //last!
    Invalid
};

} // namespace abm
} // namespace mio

#endif
