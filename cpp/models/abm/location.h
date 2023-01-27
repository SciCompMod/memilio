/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#ifndef EPI_ABM_LOCATION_H
#define EPI_ABM_LOCATION_H

#include "abm/person.h"
#include "abm/mask_type.h"
#include "abm/parameters.h"
#include "abm/location_type.h"
#include "abm/infection_state.h"
#include "abm/vaccine.h"

#include "memilio/math/eigen.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/time_series.h"
#include <array>
#include <random>

namespace mio
{
namespace abm
{
class Person;

/**
 * CellCapacity describes the size of a cell. 
 * It consists of a volume and a capacity in persons which is an upper bound for the number
 * of people that can be in the cell at the same time.
 */
struct CellCapacity {
    CellCapacity()
        : volume(0)
        , persons(std::numeric_limits<int>::max())
    {
    }
    uint32_t volume;
    uint32_t persons;
};

/**
 * The location can be split up into several Cell%s. This allows a finer division of the people in public transport.
 * By default, each Location consists of one Cell.
 */
struct Cell {
    std::vector<std::shared_ptr<Person>> m_persons;
    CustomIndexArray<ScalarType, VirusVariant, AgeGroup> m_cached_exposure_rate_contacts;
    CustomIndexArray<ScalarType, VirusVariant> m_cached_exposure_rate_air;
    CellCapacity m_capacity;

    explicit Cell(std::vector<std::shared_ptr<Person>> persons = {})
        : m_persons(std::move(persons))
        , m_cached_exposure_rate_contacts({{VirusVariant::Count, AgeGroup::Count}, 0.})
        , m_cached_exposure_rate_air({{VirusVariant::Count}, 0.})
        , m_capacity()
    {
    }

    /**
    * @brief Computes a relative cell size for the Cell.
    * @return The relative cell size for the Cell.
    */
    ScalarType compute_space_per_person_relative();

    /**
    * @brief Get subpopulation of a particular InfectionState in the Cell.
    * @param[in] t TimePoint of querry.
    * @param[in] state InfectionState of interest.
    * @return Amount of persons of the InfectionState in the Cell.
    */
    uint32_t get_subpopulation(TimePoint t, InfectionState state) const;

}; // namespace mio

/**
 * All locations in the simulated world where persons gather.
 */
class Location
{
public:
    /**
     * Construct a Location of a certain type.
     * @param type The type of the location.
     * @param index The index of the location.
     * @param num_cells [Default: 1] The number of cells in which the Location is divided.
     */
    Location(LocationType type, uint32_t index, uint32_t num_cells = 1);

    /**
     * get the type of this location.
     */
    LocationType get_type() const
    {
        return m_type;
    }

    /**
     *get the index of this location.
     */
    unsigned get_index() const
    {
        return m_index;
    }

    /**
     * @brief Compute the transmission factor for contact transmission of the virus in a Cell.
     * @param[in] cell_index Cell index of the Cell.
     * @param[in] virus VirusVariant of interest.
     * @param[in] age_receiver AgeGroup of the receiving Person.
     * @param[in] age_transmitter AgeGroup of the transmitting Person.
     * @returns Amount of average infections with the virus from the AgeGroup of the transmitter per day.
    */
    ScalarType transmission_contacts_per_day(uint32_t cell_index, VirusVariant virus, AgeGroup age_receiver,
                                             AgeGroup age_transmitter) const;

    /**
     * @brief Compute the transmission factor for a aerosol transmission of the virus in a Cell.
     * @param[in] cell_index Cell index of the Cell.
     * @param[in] virus VirusVariant of interest.
     * @returns Amount of average infections with the virus per day.
    */
    ScalarType transmission_air_per_day(uint32_t cell_index, VirusVariant virus) const;

    /** 
     * @brief A person interacts with the population at this Location and may become infected.
     * @param person The Person that interacts with the population.
     * @param dt Length of the current simulation time step.
     * @param global_params Global infection parameters.
     */
    void interact(Person& person, TimePoint t, TimeSpan dt, GlobalInfectionParameters& global_params) const;

    /** 
     * @brief Add a Person to the population at this Location.
     * @param person The Person arriving.
     * @param cell_idx [Default: 0] Index of the Cell the Person shall go to.
    */
    void add_person(const std::shared_ptr<Person>& person, uint32_t cell_idx = 0);

    /** 
     * @brief Remove a Person from the population of this Location.
     * @param person The Person leaving.
     */
    void remove_person(const std::shared_ptr<Person>& person);

    /** 
     * @brief Prepare the location for the next simulation step.
     * @param t Current TimePoint of the simulation.
     * @param dt The duration of the simulation step.
     */
    void begin_step(TimePoint t, TimeSpan dt);

    /**
     * @return parameters of the infection that are specific to this location
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
     * @return All cells of the location.
    */
    const std::vector<Cell>& get_cells() const
    {
        return m_cells;
    }

    /**
     * get the type of mask that is demanded when entering the location
     * @return type of the mask 
     */
    MaskType get_required_mask() const
    {
        return m_required_mask;
    }

    void set_required_mask(MaskType type)
    {
        m_required_mask = type;
    }

    /**
     * @brief Get the contact exposure rate in the cell.
     * @param[in] cell_idx CellIndex of interest.
     * @return Air exposure rate in the cell.
     */
    CustomIndexArray<ScalarType, VirusVariant, AgeGroup> get_cached_exposure_rate_contacts(uint32_t cell_idx)
    {
        return m_cells[cell_idx].m_cached_exposure_rate_contacts;
    }

    /**
     * @brief Get the air exposure rate in the cell.
     * @param[in] cell_idx CellIndex of interest.
     * @return Contact exposure rate in the cell.
     */
    CustomIndexArray<ScalarType, VirusVariant> get_cached_exposure_rate_air(uint32_t cell_idx)
    {
        return m_cells[cell_idx].m_cached_exposure_rate_air;
    }

    /**
    * @brief Set the capacity of a cell in the Location in persons and volume.
    * @param persons Maximum number of Person%s that can visit the Cell at the same time.
    * @param volume Volume of the Cell in m^3.
    */
    void set_capacity(uint32_t persons, uint32_t volume, uint32_t cell_idx = 0)
    {
        m_cells[cell_idx].m_capacity.persons = persons;
        m_cells[cell_idx].m_capacity.volume  = volume;
    }

    /**
    * @return The capacity of a Cell in persons and volume.
    */
    CellCapacity get_capacity(uint32_t cell_idx = 0)
    {
        return m_cells[cell_idx].m_capacity;
    }

    /**
    * Set the capacity adapted transmission risk flag
    * @param consider_capacity if true considers the capacity of the location for the computation of relative 
    * transmission risk
    */
    void set_capacity_adapted_transmission_risk_flag(bool consider_capacity)
    {
        m_capacity_adapted_transmission_risk = consider_capacity;
    }

    bool get_npi_active() const
    {
        return m_npi_active;
    }

    void set_npi_active(bool new_status)
    {
        m_npi_active = new_status;
    }

    /**
     * @brief Get the total number of Person%s at the Location.
     * @return Number of Person%s.
     */
    uint32_t get_number_persons();

    /**
    * @brief Get number of Person%s of a particular InfectionState for one Cell.
    * @param[in] t TimePoint of querry.
    * @param[in] state InfectionState of interest.
    * @param[in] cell_idx [Default: 0] CellIndex of the Cell.
    * @return Amount of Person%s of the InfectionState in the Cell.
    */
    uint32_t get_subpopulation_cell(TimePoint t, InfectionState state, uint32_t cell_idx = 0) const;

    /**
    * @brief Get the number of Person%s of a particular InfectionState for all Cell%s.
    * @param[in] t TimePoint of querry.
    * @param[in] state InfectionState of interest.
    * @return Amount of Person%s of the InfectionState in all Cell%s.
    */
    uint32_t get_subpopulation(TimePoint t, InfectionState state) const;

    /**
    * @brief Get all subpopulations for all Cell%s.
    * @param[in] t TimePoint of querry.
    * @return Vector of all subpopulations in all Cell%s.
    */
    Eigen::Ref<const Eigen::VectorXi> get_subpopulations(TimePoint t) const;

    /**
     * Add a timepoint to the subpopulations timeseries.
     * @param[in] t The TimePoint to be added.
    */
    void store_subpopulations(const TimePoint t);

    /**
     * @brief Initialize the history of subpopulations.
     * @param[in] t The TimePoint of initialization.
    */
    void initialize_subpopulations(TimePoint t);

private:
    LocationType m_type;
    uint32_t m_index;
    bool m_capacity_adapted_transmission_risk;
    LocalInfectionParameters m_parameters;
    TimeSeries<ScalarType> m_subpopulations{Eigen::Index(InfectionState::Count)};
    std::vector<Cell> m_cells{};
    MaskType m_required_mask;
    bool m_npi_active;
};

} // namespace abm
} // namespace mio

#endif
