/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Kilian Volmer
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

#include "memilio/io/default_serialize.h"
#include <cmath>
#include <math.h>
#include <boost/geometry/algorithms/distance.hpp>
#include "memilio/utils/logging.h"
namespace mio
{
namespace geo
{

const double earth_radius = 6371;
const double radians      = M_PI / 180.0;

class GeographicalLocation
{
public:
    double latitude;
    double longitude;

    /**
     * @brief Compare two GeographicalLocation%s.
     */
    bool operator==(const GeographicalLocation& other) const
    {
        return (latitude == other.latitude && longitude == other.longitude);
    }

    bool operator!=(const GeographicalLocation& other) const
    {
        return !(latitude == other.latitude && longitude == other.longitude);
    }

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("GeographicalLocation").add("latitude", latitude).add("longitude", longitude);
    }

    /*
    * @brief Calculate the distance between two geographical locations
    * @param other The other geographical location.
    * @return The distance between the two locations in kilometers.
    * 
    * Uses the haversine formula (https://en.wikipedia.org/wiki/Haversine_formula) to calculate the distance between 
    * two geographical locations. Uses an earth radius of 6375 km (https://en.wikipedia.org/wiki/Earth_radius).
    * The math functions use radians, whereas coordinates are given in degrees.
    */
    double distance(const GeographicalLocation& other) const
    {
        double delta_latitude  = (latitude - other.latitude) * radians;
        double delta_longitude = (longitude - other.longitude) * radians;
        double first_part      = sin(delta_latitude * 0.5) * sin(delta_latitude * 0.5);
        double second_part     = cos(latitude * radians) * cos(other.latitude * radians) * sin(delta_longitude * 0.5) *
                             sin(delta_longitude * 0.5);
        double distance = 2.0 * earth_radius * asin(sqrt(first_part + second_part));
        return distance;
    }
};

} // namespace geo

} // namespace mio