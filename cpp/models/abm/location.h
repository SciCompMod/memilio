/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, Khoa Nguyen, David Kerkmann
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
#ifndef MIO_ABM_LOCATION_H
#define MIO_ABM_LOCATION_H

#include "abm/location_id.h"
#include "abm/mask_type.h"
#include "abm/parameters.h"
#include "abm/location_type.h"

#include "memilio/io/auto_serialize.h"
#include "boost/atomic/atomic.hpp"

namespace mio
{
namespace abm
{

struct GeographicalLocation {
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

    /// This method is used by the auto-serialization feature.
    auto auto_serialize()
    {
        return make_auto_serialization("GraphicalLocation", NVP("latitude", latitude), NVP("longitude", longitude));
    }
};

struct CellIndex : public mio::Index<CellIndex> {
    CellIndex(size_t i)
        : mio::Index<CellIndex>(i)
    {
    }
};

using ContactExposureRates = CustomIndexArray<boost::atomic<ScalarType>, CellIndex, VirusVariant, AgeGroup>;
using AirExposureRates     = CustomIndexArray<boost::atomic<ScalarType>, CellIndex, VirusVariant>;

/**
 * @brief CellCapacity describes the size of a Cell. 
 * It consists of a volume and a capacity in Person%s which is an upper bound for the number
 * of people that can be in the Cell at the same time.
 */
struct CellCapacity {
    CellCapacity()
        : volume(0)
        , persons(std::numeric_limits<int>::max())
    {
    }
    uint32_t volume; ///< Volume of the Cell.
    uint32_t persons; ///< Maximal number of Person%s at the Cell.

    /// This method is used by the auto-serialization feature.
    auto auto_serialize()
    {
        return make_auto_serialization("CellCapacity", NVP("volume", volume), NVP("persons", persons));
    }
};

/**
 * @brief The Location can be split up into several Cell%s. 
 * This allows a finer division of the people at the Location.
 */
struct Cell {
    CellCapacity m_capacity;

    /**
    * @brief Computes a relative cell size for the Cell.
    * @return The relative cell size for the Cell.
    */
    ScalarType compute_space_per_person_relative() const;

    /// This method is used by the auto-serialization feature.
    auto auto_serialize()
    {
        return make_auto_serialization("Cell", NVP("capacity", m_capacity));
    }
}; // namespace mio

/**
 * @brief All Location%s in the simulated World where Person%s gather.
 */
class Location
{
public:
    /**
     * @brief Construct a Location with provided parameters. 
     * @param[in] loc_type The #LocationType.
     * @param[in] loc_id The index of the Location in the World.
     * @param[in] num_agegroups [Default: 1] The number of age groups in the model.
     * @param[in] num_cells [Default: 1] The number of Cell%s in which the Location is divided.
     */
    explicit Location(LocationType loc_type, LocationId loc_id, size_t num_agegroups = 1, uint32_t num_cells = 1);

    /**
     * @brief Construct a copy of a Location with a new ID.
     * @param[in] other The Location to copy from.
     * @param[in] id The ID for the new Location.
     */
    explicit Location(const Location& other, LocationId id)
        : Location(other)
    {
        m_id = id;
    }

    /**
     * @brief Compare two Location%s.
     */
    bool operator==(const Location& other) const
    {
        return (m_id == other.m_id);
    }

    bool operator!=(const Location& other) const
    {
        return !(*this == other);
    }

    /**
     * @brief Get the type of this Location.
     * @return The #LocationType of the Location.
     */
    LocationType get_type() const
    {
        return m_type;
    }

    /**
     * @brief Get the location's identifier in a World.
     * @return The location's LocationId by value.
     */
    LocationId get_id() const
    {
        return m_id;
    }

    /**
     * @brief Get the Location specific Infection parameters.
     * @return Parameters of the Infection that are specific to this Location.
     * @{
     */
    LocalInfectionParameters& get_infection_parameters()
    {
        return m_parameters;
    }

    const LocalInfectionParameters& get_infection_parameters() const
    {
        return m_parameters;
    }
    /** @} */

    /**
     * @brief Get the Cell%s of this Location.
     * @return The vector of all Cell%s of the Location.
     */
    const std::vector<Cell>& get_cells() const
    {
        return m_cells;
    }

    /**
     * @brief Get the type of Mask that is demanded when entering this Location.
     * @return Least secure MaskType that is demanded when entering this Location.
     */
    MaskType get_required_mask() const
    {
        return m_required_mask;
    }

    /**
     * @brief Set the required type of mask for entering this Location.
     * @param[in] type The Least secure MaskType that is demanded when entering this Location.
     */
    void set_required_mask(MaskType type)
    {
        m_required_mask = type;
    }

    /**
     * @brief Set the CellCapacity of a Cell in the Location in persons and volume.
     * @param[in] persons Maximum number of Person%s that can visit the Cell at the same time.
     * @param[in] volume Volume of the Cell in m^3.
     * @param[in] cell_idx Index of the Cell.
     */
    void set_capacity(uint32_t persons, uint32_t volume, uint32_t cell_idx = 0)
    {
        assert(cell_idx < m_cells.size() && "Given cell index is too large.");
        m_cells[cell_idx].m_capacity.persons = persons;
        m_cells[cell_idx].m_capacity.volume  = volume;
    }

    /**
     * @brief Get the capacity of a specific Cell in persons and volume.
     * @param[in] cell_idx The index of the Cell.
     * @return The CellCapacity of the Cell.
     */
    CellCapacity get_capacity(uint32_t cell_idx = 0) const
    {
        assert(cell_idx < m_cells.size() && "Given cell index is too large.");
        return m_cells[cell_idx].m_capacity;
    }

    /**
     * @brief Get the information whether NPIs are active at this Location.
     * If true requires e.g. Mask%s when entering a Location.
     * @return True if NPIs are active at this Location.
     */
    bool get_npi_active() const
    {
        return m_npi_active;
    }

    /**
     * @brief Activate or deactivate NPIs at this Location.
     * @param[in] new_status Status of NPIs.
     */
    void set_npi_active(bool new_status)
    {
        m_npi_active = new_status;
    }

    /**
     * @brief Get the geographical location of the Location.
     * @return The geographical location of the Location.
     */
    GeographicalLocation get_geographical_location() const
    {
        return m_geographical_location;
    }

    /**
     * @brief Set the geographical location of the Location.
     * @param[in] location The geographical location of the Location.
     */
    void set_geographical_location(GeographicalLocation location)
    {
        m_geographical_location = location;
    }

    /// This method is used by the auto-serialization feature.
    auto auto_serialize()
    {
        return make_auto_serialization("Location", NVP("id", m_id), NVP("parameters", m_parameters),
                                       NVP("cells", m_cells), NVP("required_mask", m_required_mask),
                                       NVP("npi_active", m_npi_active),
                                       NVP("geographical_location", m_geographical_location));
    }

private:
    friend AutoSerializableFactory<Location>;
    Location() = default;

    LocationType m_type; ///< Type of the Location.
    LocationId m_id; ///< Unique identifier for the Location in the World owning it.
    LocalInfectionParameters m_parameters; ///< Infection parameters for the Location.
    std::vector<Cell> m_cells{}; ///< A vector of all Cell%s that the Location is divided in.
    MaskType m_required_mask; ///< Least secure type of Mask that is needed to enter the Location.
    bool m_npi_active; ///< If true requires e.g. Mask%s to enter the Location.
    GeographicalLocation m_geographical_location; ///< Geographical location (longitude and latitude) of the Location.
};

} // namespace abm
} // namespace mio

#endif
