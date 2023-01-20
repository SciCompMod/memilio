/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Elisabeth Kluth
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

#include "abm/mask_type.h"
#include "abm/parameters.h"
#include "abm/state.h"
#include "abm/location_type.h"

#include "memilio/math/eigen.h"
#include "memilio/utils/custom_index_array.h"
#include <array>
#include <random>

namespace mio
{
namespace abm
{
class Person;

/**
 * @brief LocationCapacity describes the size of a Location. 
 * It consists of a volume and a capacity in Person%s which is an upper bound for the number
 * of Person%s that can be at the Location at the same time.
 */
struct LocationCapacity {
    LocationCapacity()
        : volume(0)
        , persons(std::numeric_limits<int>::max())
    {
    }
    int volume; ///< Capacity in volume.
    int persons; ///< Capactiy in Person%s.
};

/**
 * @brief LocationId identifies a Location uniquely.
 * It consists of the LocationType of the Location and an Index. The index corresponds to the index into the structure
 * m_locations from World, where all Location%s are saved.
 */
struct LocationId {
    uint32_t index; ///< Unique index of the Location.
    LocationType type; ///< Type of the Location.

    bool operator==(const LocationId& rhs) const
    {
        return (index == rhs.index && type == rhs.type);
    }

    bool operator!=(const LocationId& rhs) const
    {
        return !(index == rhs.index && type == rhs.type);
    }
};

/**
 * @brief A finer division of the Location.
 * The Location can be split up into several Cell%s. This allows a finer division of the people at the Location.
 */
struct Cell {
    uint32_t num_people; ///< Number of Person%s in the Cell.
    uint32_t num_carriers; ///< Number of pre- and asymptomatic Person%s in the Cell.
    uint32_t num_infected; ///< Number of symptomatic Person%s in the Cell.
    /** @todo describe exposure rate*/
    CustomIndexArray<double, AgeGroup, VaccinationState> cached_exposure_rate; ///<

    Cell()
        : num_people(0)
        , num_carriers(0)
        , num_infected(0)
        , cached_exposure_rate({{AgeGroup::Count, VaccinationState::Count}, 0.})
    {
    }

    Cell(uint32_t num_p, uint32_t num_c, uint32_t num_i,
         CustomIndexArray<double, AgeGroup, VaccinationState> cached_exposure_rate_new)
        : num_people(num_p)
        , num_carriers(num_c)
        , num_infected(num_i)
        , cached_exposure_rate(cached_exposure_rate_new)
    {
    }

}; // namespace mio

/**
 * All Location%s in the simulated World where Person%s gather.
 */
class Location
{
public:
    /**
     * @brief Construct a Location of a certain type.
     * @param[in] type The type of the Location.
     * @param[in] index The index of the Location.
     * @param[in] num_cells The number of Cell%s in which the Location is divided.
     */
    Location(LocationType type, uint32_t index, uint32_t num_cells = 0);

    /**
     * @brief Get the type of this Location.
     * @return The type of this Location.
     */
    LocationType get_type() const
    {
        return m_type;
    }

    /**
     * @brief Get the index of this Location.
     * @return The index of this Location.
     */
    unsigned get_index() const
    {
        return m_index;
    }

    /** 
     * @brief A Person interacts with the population at this Location which may change the #InfectionState.
     * @param[in] person The Person that interacts with the population.
     * @param[in] dt Length of the current Simulation time step.
     * @param[in] global_params Infection parameters that are the same everywhere in the World.
     * @return New #InfectionState of the Person.
     */
    InfectionState interact(const Person& person, TimeSpan dt, const GlobalInfectionParameters& global_params) const;

    /** 
     * @brief Add a Person to the population at this Location.
     * @param[in] person The Person arriving.
    */
    void add_person(const Person& person);

    /** 
     * @brief Remove a Person from the population of this Location.
     * @param[in] person The Person leaving.
     */
    void remove_person(const Person& person);

    /** 
     * @brief Notification that one Person in this Location changed the #InfectionState.
     * @param[in] person The Person that changed the #InfectionState.
     * @param[in] old_state The previous #InfectionState of the Person.
     */
    void changed_state(const Person& person, InfectionState old_infection_state);

    /** 
     * @brief Prepare the Location for the next Simulation step.
     * @param[in] dt The duration of the Simulation step.
     * @param[in] global_params Infection parameters that are the same everywhere in the World.
     */
    void begin_step(TimeSpan dt, const GlobalInfectionParameters& global_params);

    /** 
     * @brief Get the number of Person%s at this Location in one #InfectionState.
     * @return The number of Person%s at this Location that are in the specified #InfectionState.
     */
    int get_subpopulation(InfectionState s) const;

    /** 
     * @brief Get the number of Person%s at this Location for all #InfectionState%s.
     * The vector is indexed by #InfectionState.
     * @return The number of Person%s in all #InfectionState%s.
     * */
    Eigen::Ref<const Eigen::VectorXi> get_subpopulations() const;

    /**
     * @brief Get the Location specific infection parameters.
     * @return parameters of the infection that are specific to this Location.
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
     * @return A vector with the Cell%s.
     */
    const std::vector<Cell>& get_cells() const
    {
        return m_cells;
    }

    /**
     * @brief Get the type of Mask that is demanded when entering this Location.
     * @return The type of the Mask.
     */
    MaskType get_required_mask() const
    {
        return m_required_mask;
    }

    /**
     * @brief Set the required type of Mask for entering this Location.
     * @param[in] type The type of the Mask.
     */
    void set_required_mask(MaskType type)
    {
        m_required_mask = type;
    }

    /**
     * @brief Get the number of Person%s at this Location.
     * @return The number of Person%s.
     */
    int get_population()
    {
        return m_num_persons;
    }

    /**
     * @brief Get the exposure rate of the Location.
     * @return The exposure rate for every #AgeGroup and #VaccinationState.
     */
    CustomIndexArray<double, AgeGroup, VaccinationState> get_cached_exposure_rate()
    {
        return m_cached_exposure_rate;
    }

    /**
     * @brief Set the capacity of the Location in Person and volume.
     * @param[in] persons The maximum number of Person%s that can visit the Location at the same time.
     * @param[in] volume The volume of the Location in m^3.
     */
    void set_capacity(int persons, int volume)
    {
        m_capacity.persons = persons;
        m_capacity.volume  = volume;
    }

    /**
     * @brief Get the capacity of the Location in Person and volume.
     * @return The capacity in Person and volume.
     */
    LocationCapacity get_capacity()
    {
        return m_capacity;
    }

    /**
     * @brief Computes a relative transmission risk factor for the Location.
     * @return The relative risk factor for the Location.
     */
    double compute_relative_transmission_risk();

    /**
     * @brief Set the capacity adapted transmission risk flag.
     * @param[in] consider_capacity If true considers the capacity of the location for the computation of relative 
     * transmission risk.
     */
    void set_capacity_adapted_transmission_risk_flag(bool consider_capacity)
    {
        m_capacity_adapted_transmission_risk = consider_capacity;
    }

    /**
     * @brief Get the information if NPIs are active at this Location.
     * If true requires e.g. Mask%s when entering a Location.
     * @return True if NPIs are active.
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

private:
    /**
     * @brief Update the number of Person%s in an #InfectionState at this Location.
     * @param[in] s The #InfectionState of which the number of Person%s has changed.
     * @param[in] delta The change in number of Person%s.
     */
    void change_subpopulation(InfectionState s, int delta);

private:
    LocationType m_type;
    uint32_t m_index;
    int m_num_persons = 0;
    LocationCapacity m_capacity;
    bool m_capacity_adapted_transmission_risk;
    std::array<int, size_t(InfectionState::Count)> m_subpopulations;
    LocalInfectionParameters m_parameters;
    CustomIndexArray<double, AgeGroup, VaccinationState> m_cached_exposure_rate;
    std::vector<Cell> m_cells;
    MaskType m_required_mask;
    bool m_npi_active;
};

} // namespace abm
} // namespace mio

#endif
