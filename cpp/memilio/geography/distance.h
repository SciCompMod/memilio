/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Kilian Volmer, Rene Schmiedling
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
#ifndef MIO_DISTANCE_H
#define MIO_DISTANCE_H

#include "memilio/config.h"
#include "memilio/io/default_serialize.h"

namespace mio
{

namespace geo
{

/**
 * @brief Represents a point in time.
 * Seconds from an unspecified monday at 00:00 (epoch). 
 */
class Distance
{
public:
    /**
     * @brief Default ctor, unitinialized.
     */
    Distance() = default;
    /**
     * @brief Creates a Distance from a specified number of seconds.
     * @param[in] seconds The number of seconds after the epoch.
     */
    explicit Distance(ScalarType meters)
        : m_meters(meters)
    {
    }

    /**
     * @brief Distance in kilometers.
     */
    ScalarType kilometers() const
    {
        return ScalarType(m_meters) / 1000;
    }

    /**
     * @brief Distance in meters.
     */
    ScalarType meters() const
    {
        return m_meters;
    }

    /**
     * @name Comparison operators.
     * @{
     */
    bool operator==(const Distance& other) const
    {
        return m_meters == other.m_meters;
    }
    bool operator!=(const Distance& other) const
    {
        return !(*this == other);
    }
    bool operator<(const Distance& other) const
    {
        return m_meters < other.m_meters;
    }
    bool operator<=(const Distance& other) const
    {
        return m_meters <= other.m_meters;
    }
    bool operator>(const Distance& other) const
    {
        return m_meters > other.m_meters;
    }
    bool operator>=(const Distance& other) const
    {
        return m_meters >= other.m_meters;
    }
    /**@}*/

    /**
     * @brief Add or subtract a Distance.
     * @{
     */
    Distance operator+(const Distance& m) const
    {
        return Distance{m_meters + m.meters()};
    }
    Distance& operator+=(const Distance& m)
    {
        m_meters += m.meters();
        return *this;
    }
    Distance operator-(const Distance& m) const
    {
        return Distance{m_meters - m.meters()};
    }
    Distance& operator-=(const Distance& m)
    {
        m_meters -= m.meters();
        return *this;
    }
    /**@}*/

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("Distance").add("meters", m_meters);
    }

private:
    int m_meters; ///< The distance in meters.
};

/**
 * @brief Create a Distance of a specified number of meters.
 * @param[in] meters Number of meters in the Distance.
 */
inline Distance meters(ScalarType meters)
{
    return Distance(meters);
}

/**
 * @brief Create a Distance of a specified number of kilometers.
 * @param[in] kilometers Number of kilometers in the Distance.
 */
inline Distance kilometers(ScalarType kilometers)
{
    return Distance(kilometers * 1000);
}

} // namespace geo
} // namespace mio

#endif
