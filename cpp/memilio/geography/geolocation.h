/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Kilian Volmer, Sascha Korf, Carlotta Gerstein, Daniel Abele, Elisabeth Kluth, Khoa Nguyen, David Kerkmann
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
#ifndef MIO_GEOGRAPHY_LOCATIONS_H
#define MIO_GEOGRAPHY_LOCATIONS_H

#include "memilio/io/default_serialize.h"
#include <cassert>
#include <cmath>
#include <numbers>
namespace mio
{
namespace geo
{
/**
 * @brief Class representing a geographical location on the Earth's surface.
 *
 * Provides latitude and longitude in degrees and a method to calculate the distance to another location.
 */
class GeographicalLocation
{

public:
    /**
     * @brief Construct a new Geographical Location object.
     * 
     * @param lat Latitude in degrees.
     * @param lon Longitude in degrees.
     */
    GeographicalLocation(double lat, double lon)
        : m_latitude(lat)
        , m_longitude(lon)
    {
        check_input();
    }
    /**
     * @brief Structure representing geographical coordinates for constructing a GeographicalLocation object.
     * 
     * @param latitude Latitude in degrees.
     * @param longitude Longitude in degrees.
     */
    /**
     * @brief Construct a new Geographical Location object.
     * 
     * @param coordinates Pair of latitude and longitude in degrees as doubles.
     */
    GeographicalLocation(std::pair<double, double> coordinates)
        : m_latitude(coordinates.first)
        , m_longitude(coordinates.second)
    {
        check_input();
    }

    /**
     * @brief Construct a new Geographical Location object using the default constructor to keep compatibility with ABM code.
     * 
     */
    GeographicalLocation() = default;

    /**
     * @brief Compare two GeographicalLocation%s for equality
     */
    bool operator==(const GeographicalLocation& other) const
    {
        return (m_latitude == other.m_latitude && m_longitude == other.m_longitude);
    }

    /**
     * @brief Compare two GeographicalLocation%s for inequality
     */
    bool operator!=(const GeographicalLocation& other) const
    {
        return !(m_latitude == other.m_latitude && m_longitude == other.m_longitude);
    }

    /**
     * @brief Default serialize the GeographicalLocation object.
     * 
     * This method is used by the default serialization feature.
     */
    auto default_serialize()
    {
        return Members("GeographicalLocation").add("latitude", m_latitude).add("longitude", m_longitude);
    }

    void check_input() const
    {
        assert(m_latitude <= 90. && m_latitude >= -90.);
        assert(m_longitude <= 180. && m_longitude >= -180.);
    }

    /*
    * @brief Calculate the distance between two GeographicalLocation%s.
    * @param other The other GeographicalLocation.
    * @return The distance between the two locations in kilometers.
    * 
    * Uses the haversine formula (https://en.wikipedia.org/wiki/Haversine_formula) to calculate the distance between 
    * two geographical locations. Uses an earth radius of 6371 km (https://en.wikipedia.org/wiki/Earth_radius).
    * The math functions use radians, whereas coordinates are given in degrees.
    */
    double distance(const GeographicalLocation& other) const
    {
        const double delta_latitude           = (m_latitude - other.m_latitude) * radians;
        const double delta_longitude          = (m_longitude - other.m_longitude) * radians;
        const double sin_delta_latitude_half  = sin(delta_latitude * 0.5);
        const double sin_delta_longitude_half = sin(delta_longitude * 0.5);
        const double first_part               = sin_delta_latitude_half * sin_delta_latitude_half;
        const double second_part              = cos(m_latitude * radians) * cos(other.m_latitude * radians) *
                                   sin_delta_longitude_half * sin_delta_longitude_half;
        const double distance = 2.0 * earth_radius * asin(sqrt(first_part + second_part));
        return distance;
    }

    /**
     * @brief Get the latitude object
     * 
     * @return double latitude in degrees
     */
    double get_latitude() const
    {
        return m_latitude;
    }

    /**
     * @brief Get the longitude object
     * 
     * @return double longitude in degrees
     */
    double get_longitude() const
    {
        return m_longitude;
    }

    /**
     * @brief Set the latitude object
     * 
     * @param lat Latitude in degrees.
     */
    void set_latitude(double lat)
    {
        m_latitude = lat;
    }

    /**
     * @brief Set the longitude object
     * 
     * @param lon Longitude in degrees.
     */
    void set_longitude(double lon)
    {
        m_longitude = lon;
    }

private:
    double m_latitude;
    double m_longitude;
    constexpr static double earth_radius = 6371;
    constexpr static double radians      = std::numbers::pi / 180.0;
};

} // namespace geo

} // namespace mio

#endif // MIO_GEOGRAPHY_LOCATIONS_H
