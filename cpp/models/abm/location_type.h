/* 
* Copyright (C) 2020-2024 MEmilio
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
#ifndef EPI_ABM_LOCATION_TYPE_H
#define EPI_ABM_LOCATION_TYPE_H

#include <cstdint>
#include <limits>
#include <functional>

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
    Cemetery, // Location for all the dead persons. It is created once for the World.
    EventPanvadere,

    Count //last!
};

static constexpr uint32_t INVALID_LOCATION_INDEX = std::numeric_limits<uint32_t>::max();

/**
 * LocationId identifies a Location uniquely. It consists of the LocationType of the Location and an Index.
 * The index corresponds to the index into the structure m_locations from world, where all Locations are saved.
 */
struct LocationId {
    uint32_t index;
    LocationType type;

    bool operator==(const LocationId& rhs) const
    {
        return (index == rhs.index && type == rhs.type);
    }

    bool operator!=(const LocationId& rhs) const
    {
        return !(index == rhs.index && type == rhs.type);
    }

    bool operator<(const LocationId& rhs) const
    {
        if (type == rhs.type) {
            return index < rhs.index;
        }
        return (type < rhs.type);
    }
};

struct GeographicalLocation {
    double latitude;
    double longitude;

    /**
     * @brief Compare two Location%s.
     */
    bool operator==(const GeographicalLocation& other) const
    {
        return (latitude == other.latitude && longitude == other.longitude);
    }

    bool operator!=(const GeographicalLocation& other) const
    {
        return !(latitude == other.latitude && longitude == other.longitude);
    }
};

} // namespace abm
} // namespace mio

template <>
struct std::hash<mio::abm::LocationId> {
    std::size_t operator()(const mio::abm::LocationId& loc_id) const
    {
        return (std::hash<uint32_t>()(loc_id.index)) ^ (std::hash<uint32_t>()(static_cast<uint32_t>(loc_id.type)));
    }
};

#endif
