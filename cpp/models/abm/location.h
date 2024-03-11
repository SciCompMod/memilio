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

#include "abm/mask_type.h"
#include "abm/parameters.h"
#include "abm/location_type.h"

#include "boost/atomic/atomic.hpp"

namespace mio
{
namespace abm
{

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

    /**
    * @brief Get subpopulation of a particular #InfectionState in the Cell.
    * @param[in] t TimePoint of querry.
    * @param[in] state #InfectionState of interest.
    * @return Amount of Person%s of the #InfectionState in the Cell.
    */
    // size_t get_subpopulation(TimePoint t, InfectionState state) const;

}; // namespace mio

/**
 * @brief All Location%s in the simulated World where Person%s gather.
 */
class Location
{
public:
    /**
     * @brief Construct a Location of a certain LocationId.
     * @param[in] loc_id The #LocationId.
     * @param[in] num_agegroups [Default: 1] The number of age groups in the model.
     * @param[in] num_cells [Default: 1] The number of Cell%s in which the Location is divided.
     */
    explicit Location(LocationId loc_id, size_t num_agegroups = 1, uint32_t num_cells = 1);

    /**
     * @brief Construct a Location with provided parameters. 
     * @param[in] loc_type The #LocationType.
     * @param[in] index The index of the Location.
     * @param[in] num_agegroups [Default: 1] The number of age groups in the model.
     * @param[in] num_cells [Default: 1] The number of Cell%s in which the Location is divided.
     */
    explicit Location(LocationType loc_type, uint32_t loc_index, size_t num_agegroups = 1, uint32_t num_cells = 1)
        : Location(LocationId{loc_index, loc_type}, num_agegroups, num_cells)
    {
    }

    /**
     * @brief Construct a copy of a Location with a new ID.
     * @param[in] other A Location.
     * @param[in] id The ID for the new Location.
     */
    explicit Location(const Location& other, LocationId id)
        : Location(other)
    {
        m_id = id;
    }

    /**
     * @brief Return a copy of this #Location object with an empty m_persons.
     * @param[in] num_agegroups The number of age groups in the model.
     */
    Location copy() const;

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
        return m_id.type;
    }

    /**
     * @brief Get the index of this Location.
     * @return The index of the Location.
     */
    unsigned get_index() const
    {
        return m_id.index;
    }

    /**
     * @brief Compute the transmission factor for contact transmission of the virus in a Cell.
     * @param[in] cell_index Cell index of the Cell.
     * @param[in] virus VirusVariant of interest.
     * @param[in] age_receiver AgeGroup of the receiving Person.
     * @param[in] age_transmitter AgeGroup of the transmitting Person.
     * @param[in] num_agegroups The number of age groups in the model.
     * @return Amount of average Infection%s with the virus from the AgeGroup of the transmitter per day.
    */
    // ScalarType transmission_contacts_per_day(uint32_t cell_index, VirusVariant virus, AgeGroup age_receiver,
    //                                          size_t num_agegroups) const;

    /**
     * @brief Compute the transmission factor for a aerosol transmission of the virus in a Cell.
     * @param[in] cell_index Cell index of the Cell.
     * @param[in] virus VirusVariant of interest.
     * @param[in] global_params The Parameters set of the World. 
     * @return Amount of average Infection%s with the virus per day.
    */
    // ScalarType transmission_air_per_day(uint32_t cell_index, VirusVariant virus, const Parameters& global_params) const;

    /** 
     * @brief Prepare the Location for the next Simulation step.
     * @param[in] t Current TimePoint of the Simulation.
     * @param[in] dt The duration of the Simulation step.
     * @param[in] num_agegroups The number of age groups in the model.
     */
    // void cache_exposure_rates(TimePoint t, TimeSpan dt, size_t num_agegroups);

    /**
     * @brief Get the Location specific Infection parameters.
     * @return Parameters of the Infection that are specific to this Location.
     */
    LocalInfectionParameters& get_infection_parameters()
    {
        return m_parameters;
    }

    const LocalInfectionParameters& get_infection_parameters() const
    {
        return m_parameters;
    }

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
        assert(cell_idx < m_cells.size());
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
        assert(cell_idx < m_cells.size());
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
     * serialize this. 
     * @see mio::serialize
     */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Location");
        obj.add_element("index", m_id.index);
        obj.add_element("type", m_id.type);
    }

    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Location> deserialize(IOContext& io)
    {
        auto obj   = io.expect_object("Location");
        auto index = obj.expect_element("index", Tag<uint32_t>{});
        auto type  = obj.expect_element("type", Tag<uint32_t>{});
        return apply(
            io,
            [](auto&& index_, auto&& type_) {
                return Location{LocationId{index_, LocationType(type_)}};
            },
            index, type);
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

    // return id by value. used to identify a location in a World
    LocationId get_id() const
    {
        return m_id;
    }

private:
    LocationId m_id; ///< Id of the Location including type and index.
    LocalInfectionParameters m_parameters; ///< Infection parameters for the Location.
    std::vector<Cell> m_cells{}; ///< A vector of all Cell%s that the Location is divided in.
    MaskType m_required_mask; ///< Least secure type of Mask that is needed to enter the Location.
    bool m_npi_active; ///< If true requires e.g. Mask%s to enter the Location.
    GeographicalLocation m_geographical_location; ///< Geographical location (longitude and latitude) of the Location.
};

} // namespace abm
} // namespace mio

#endif
