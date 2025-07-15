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
#ifndef EPI_ABM_LOCATION_H
#define EPI_ABM_LOCATION_H

#include "abm/person.h"
#include "abm/mask_type.h"
#include "abm/parameters.h"
#include "abm/location_type.h"
#include "abm/infection_state.h"
#include "abm/vaccine.h"
#include "abm/time.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/memory.h"
#include "memilio/utils/random_number_generator.h"
#include <array>
#include <random>
#include <mutex>

namespace mio
{
namespace abm
{
class Person;

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
    std::vector<observer_ptr<Person>> m_persons;
    CustomIndexArray<ScalarType, VirusVariant, AgeGroup> m_cached_exposure_rate_contacts;
    CustomIndexArray<ScalarType, VirusVariant> m_cached_exposure_rate_air;
    CellCapacity m_capacity;

    Cell(size_t num_agegroups, std::vector<observer_ptr<Person>> persons = {})
        : m_persons(std::move(persons))
        , m_cached_exposure_rate_contacts({{VirusVariant::Count, AgeGroup(num_agegroups)}, 0.})
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
    * @brief Get subpopulation of a particular #InfectionState in the Cell.
    * @param[in] t TimePoint of querry.
    * @param[in] state #InfectionState of interest.
    * @return Amount of Person%s of the #InfectionState in the Cell.
    */
    size_t get_subpopulation(TimePoint t, InfectionState state) const;

}; // namespace mio

/**
 * @brief All Location%s in the simulated World where Person%s gather.
 */
using HourlyContactMatrix = std::array<Eigen::MatrixXd, 24>;
class Location
{
public:
    /**
     * @brief Construct a Location of a certain LocationId.
     * @param[in] loc_id The #LocationId.
     * @param[in] num_agegroups [Default: 1] The number of age groups in the model.
     * @param[in] num_cells [Default: 1] The number of Cell%s in which the Location is divided.
     */
    Location(LocationId loc_id, size_t num_agegroups = 1, uint32_t num_cells = 1);

    /**
     * @brief Construct a Location with provided parameters. 
     * @param[in] loc_type The #LocationType.
     * @param[in] index The index of the Location.
     * @param[in] num_agegroups [Default: 1] The number of age groups in the model.
     * @param[in] num_cells [Default: 1] The number of Cell%s in which the Location is divided.
     */
    Location(LocationType loc_type, uint32_t loc_index, size_t num_agegroups = 1, uint32_t num_cells = 1)
        : Location(LocationId{loc_index, loc_type}, num_agegroups, num_cells)
    {
    }

    /**
     * @brief Return a copy of the current Location object.
     * @param[in] other The original #Location.
     */
    Location(const Location& other)
        : m_id(other.m_id)
        , m_capacity_adapted_transmission_risk(other.m_capacity_adapted_transmission_risk)
        , m_parameters(other.m_parameters)
        , m_persons(other.m_persons)
        , m_cells(other.m_cells)
        , m_required_mask(other.m_required_mask)
        , m_npi_active(other.m_npi_active)
        , m_geographical_location(other.m_geographical_location)
    {
    }

    bool location_contaminated = false;

    /**
     * @brief Return a copy of this #Location object with an empty m_persons.
     * @param[in] num_agegroups The number of age groups in the model.
     */
    Location copy_location_without_persons(size_t num_agegroups);

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
    ScalarType transmission_contacts_per_day(uint32_t cell_index, VirusVariant virus, AgeGroup age_receiver,
                                             size_t num_agegroups) const;

    /**
     * @brief Compute the transmission factor for a aerosol transmission of the virus in a Cell.
     * @param[in] cell_index Cell index of the Cell.
     * @param[in] virus VirusVariant of interest.
     * @param[in] global_params The Parameters set of the World. 
     * @return Amount of average Infection%s with the virus per day.
    */
    ScalarType transmission_air_per_day(uint32_t cell_index, VirusVariant virus, const Parameters& global_params) const;

    /** 
     * @brief A Person interacts with the population at this Location and may become infected. But for the micro matrices.
     * @param[in, out] rng Person::RandomNumberGenerator for this Person.
     * @param[in, out] person The Person that interacts with the population.
     * @param[in] dt Length of the current Simulation time step.
     * @param[in] params Parameters of the Model.
     */
    void interact_micro(Person::RandomNumberGenerator& rng, Person& person, TimePoint t, TimeSpan dt,
                        const Parameters& params) const;

    /** 
     * @brief A Person interacts with the population at this Location and may become infected.
     * @param[in, out] rng Person::RandomNumberGenerator for this Person.
     * @param[in, out] person The Person that interacts with the population.
     * @param[in] dt Length of the current Simulation time step.
     * @param[in] params Parameters of the Model.
     */
    void interact(Person::RandomNumberGenerator& rng, Person& person, TimePoint t, TimeSpan dt,
                  const Parameters& params);

    /** 
     * @brief Add a Person to the population at this Location.
     * @param[in] person The Person arriving.
     * @param[in] cell_idx [Default: 0] Index of the Cell the Person shall go to.
    */
    void add_person(Person& person, std::vector<uint32_t> cells = {0});

    /** 
     * @brief Remove a Person from the population of this Location.
     * @param[in] person The Person leaving.
     */
    void remove_person(Person& person);

    /**
     * @brief Adjust local contact rates based on maximum contacts.
     * @param[in] num_agegroups The number of age groups in the model.
    */
    void adjust_contact_rates(size_t num_agegroups);

    /** 
     * @brief Prepare the Location for the next Simulation step.
     * @param[in] t Current TimePoint of the Simulation.
     * @param[in] dt The duration of the Simulation step.
     * @param[in] num_agegroups The number of age groups in the model.
     * @param[in] params Parameters of the Model.
     */
    void cache_exposure_rates(TimePoint t, TimeSpan dt, size_t num_agegroups);

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
     * @brief add a damping (chance to enter the location) to the location
     * @return The LocationId of the Location.
     */
    void add_damping(TimePoint t, double p);

    bool entry_allowed_dampings(Person::RandomNumberGenerator& rng, const mio::abm::TimePoint t);

    /**
     * @brief Set the required type of mask for entering this Location.
     * @param[in] type The Least secure MaskType that is demanded when entering this Location.
     */
    void set_required_mask(MaskType type)
    {
        m_required_mask = type;
    }

    /**
     * @brief Get the contact exposure rate in the Cell.
     * @param[in] cell_idx Cell index of interest.
     * @return Air exposure rate in the Cell.
     */
    CustomIndexArray<ScalarType, VirusVariant, AgeGroup> get_cached_exposure_rate_contacts(uint32_t cell_idx) const
    {
        return m_cells[cell_idx].m_cached_exposure_rate_contacts;
    }

    /**
     * @brief Get the air exposure rate in the Cell.
     * @param[in] cell_idx Cell index of interest.
     * @return Contact exposure rate in the cell.
     */
    CustomIndexArray<ScalarType, VirusVariant> get_cached_exposure_rate_air(uint32_t cell_idx) const
    {
        return m_cells[cell_idx].m_cached_exposure_rate_air;
    }

    /**
     * @brief Set the CellCapacity of a Cell in the Location in persons and volume.
     * @param[in] persons Maximum number of Person%s that can visit the Cell at the same time.
     * @param[in] volume Volume of the Cell in m^3.
     * @param[in] cell_idx Index of the Cell.
     */
    void set_capacity(uint32_t persons, uint32_t volume, uint32_t cell_idx = 0)
    {
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
        return m_cells[cell_idx].m_capacity;
    }

    /**
     * @brief Set the capacity adapted transmission risk flag.
     * @param[in] consider_capacity If true considers the capacity of the Cell%s of this Location for the computation of 
     * relative transmission risk.
     */
    void set_capacity_adapted_transmission_risk_flag(bool consider_capacity)
    {
        m_capacity_adapted_transmission_risk = consider_capacity;
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
     * @brief Get the total number of Person%s at the Location.
     * @return Number of Person%s.
     */
    size_t get_number_persons() const;

    /**
     * @brief Get the number of Person%s of a particular #InfectionState for all Cell%s.
     * @param[in] t TimePoint of querry.
     * @param[in] state #InfectionState of interest.
     * @return Amount of Person%s of the #InfectionState in all Cell%s.
     */
    size_t get_subpopulation(TimePoint t, InfectionState state) const;

    /**
     * @brief Get the number of Person%s of a particular #InfectionState for all Cell%s.
     * @param[in] t TimePoint of querry.
     * @param[in] state #InfectionState of interest.
     * @param[in] age_group AgeGroup of interest.
     * @return Amount of Person%s of the #InfectionState in all Cell%s.
     */
    size_t get_subpopulation_per_age_group(TimePoint t, InfectionState state, AgeGroup age_group) const;

    /**
     * @brief Get the geographical location of the Location.
     * @return The geographical location of the Location.
     */
    GeographicalLocation get_geographical_location() const
    {
        return m_geographical_location;
    }
    void clear_infected_persons()
    {
        m_track_infected_persons.clear();
    }

    std::vector<uint32_t> get_infected_persons() const
    {
        return m_track_infected_persons;
    }

    /**
     * @brief Set the geographical location of the Location.
     * @param[in] location The geographical location of the Location.
     */
    void set_geographical_location(GeographicalLocation location)
    {
        m_geographical_location = location;
    }

    void assign_contact_matrices(HourlyContactMatrix contact_matrices, std::vector<uint32_t> assigned_persons,
                                 bool dynamic_assignment)
    {
        m_hourly_contact_matrices = contact_matrices;
        m_assigned_persons        = assigned_persons;
        m_dynamic_assignment      = dynamic_assignment;
    }

    bool is_contact_matrix_assignment_dynamic() const
    {
        return m_dynamic_assignment;
    }

    const HourlyContactMatrix& get_contact_matrices() const
    {
        return m_hourly_contact_matrices;
    }

    const std::vector<uint32_t>& get_assigned_persons() const
    {
        return m_assigned_persons;
    }

    const auto& get_persons() const
    {
        return m_persons;
    }

private:
    std::mutex m_mut; ///< Mutex to protect the list of persons from concurrent modification.
    LocationId m_id; ///< Id of the Location including type and index.
    bool m_capacity_adapted_transmission_risk; /**< If true considers the LocationCapacity for the computation of the 
    transmission risk.*/
    LocalInfectionParameters m_parameters; ///< Infection parameters for the Location.
    std::vector<std::pair<mio::abm::TimePoint, double>> m_npi_damping; ///< Temporary storage for NPI dampings.
    std::vector<observer_ptr<Person>> m_persons{}; ///< A vector of all Person%s at the Location.
    std::vector<Cell> m_cells{}; ///< A vector of all Cell%s that the Location is divided in.
    MaskType m_required_mask; ///< Least secure type of Mask that is needed to enter the Location.
    bool m_npi_active; ///< If true requires e.g. Mask%s to enter the Location.
    GeographicalLocation m_geographical_location; ///< Geographical location (longitude and latitude) of the Location.
    HourlyContactMatrix m_hourly_contact_matrices; ///< Contact matrices for the Location.
    mutable std::vector<uint32_t> m_assigned_persons;
    bool m_dynamic_assignment = false;
    std::vector<uint32_t> m_track_infected_persons; ///< List of infected Person%s at the Location this timestep.
};

} // namespace abm
} // namespace mio

#endif
