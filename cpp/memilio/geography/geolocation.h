/* 
* Copyright (C) 2020-2026 MEmilio
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

#include "distance.h"
#include "memilio/io/default_serialize.h"
#include "memilio/geography/distance.h"
#include "memilio/config.h"
#include "memilio/utils/logging.h"
#include <cassert>
#include <cmath>
#include <numbers>
namespace mio
{
namespace geo
{
/**
 * @brief Class representing a geographical location in a x-y coordinate system.
 *
 * Provides latitude and longitude in meters and a method to calculate the distance to another location.
 */
class GeographicalLocation
{

public:
    /**
     * @brief Construct a new Geographical Location object.
     * 
     * @param x X coordinate in meters.
     * @param y Y coordinate in meters.
     */
    GeographicalLocation(ScalarType x, ScalarType y)
        : m_latitude(x)
        , m_longitude(y)
    {
    }
    /**
     * @brief Construct a new Geographical Location object.
     * 
     * @param coordinates Pair of latitude and longitude in degrees as ScalarTypes.
     */
    GeographicalLocation(std::pair<ScalarType, ScalarType> coordinates)
        : m_latitude(coordinates.first)
        , m_longitude(coordinates.second)
    {
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
        return !(*this == other);
    }

    /**
     * @brief Check that this location is within a (small) distance of another.
     * @param[in] other Any location.
     * @param[in] tol The Absolute tolerance for considering two locations close.
     */
    bool is_close(const GeographicalLocation& other,
                  Distance tol = Distance(10 * Limits<ScalarType>::zero_tolerance())) const
    {
        return distance(other) < tol;
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

    /*
    * @brief Calculate the distance between two GeographicalLocation%s.
    * @param other The other GeographicalLocation.
    * @return The distance between the two locations as @ref mio::geo::Distance .
    * 
    * Uses the haversine formula (https://en.wikipedia.org/wiki/Haversine_formula) to calculate the distance between 
    * two geographical locations. Uses an earth radius of 6371 km (https://en.wikipedia.org/wiki/Earth_radius).
    * The math functions use radians, whereas coordinates are given in degrees.
    */
    Distance distance(const GeographicalLocation& other) const
    {
        const auto distance_in_meters = sqrt((other.get_x() - m_latitude) * (other.get_x() - m_latitude) +
                                             (other.get_y() - m_longitude) * (other.get_y() - m_longitude));
        return meters(distance_in_meters);
    }

    /**
     * @brief Get the x coordinate
     * 
     * @return ScalarType x coordinate in meters
     */
    ScalarType get_x() const
    {
        return m_latitude;
    }

    /**
     * @brief Get the y coordinate
     * 
     * @return ScalarType y coordinate in meters
     */
    ScalarType get_y() const
    {
        return m_longitude;
    }

    /**
     * @brief Set the x coordinate object
     * 
     * @param x X coordinate in meters.
     */
    void set_x(ScalarType x)
    {
        m_latitude = x;
    }

    /**
     * @brief Set the y coordinate object
     * 
     * @param y Y coordinate in meters.
     */
    void set_y(ScalarType y)
    {
        m_longitude = y;
    }
    /**
     * @brief Get the latitude object
     * 
     * @return latitude 
     */
    auto get_latitude() const
    {
        return m_latitude;
    }

    /**
     * @brief Get the longitude object
     * 
     * @return longitude 
     */
    auto get_longitude() const
    {
        return m_longitude;
    }

private:
    ScalarType m_latitude;
    ScalarType m_longitude;
};

} // namespace geo

} // namespace mio

#endif // MIO_GEOGRAPHY_LOCATIONS_H
