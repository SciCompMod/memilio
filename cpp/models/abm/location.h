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

#include "abm/mask_type.h"
#include "abm/parameters.h"
#include "abm/location_type.h"
#include "abm/infection_state.h"
#include "abm/vaccine.h"
#include "abm/mask.h"
#include "memilio/utils/random_number_generator.h"
#include "abm/random_events.h"

#include "memilio/math/eigen.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/memory.h"
#include "abm/infection.h"
#include <array>
#include <random>
#include <numeric>

namespace mio
{
namespace abm
{
template<typename>
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
template<typename FP=double>
struct Cell {
    std::vector<observer_ptr<Person<FP>>> m_persons;
    CustomIndexArray<ScalarType, VirusVariant, AgeGroup> m_cached_exposure_rate_contacts;
    CustomIndexArray<ScalarType, VirusVariant> m_cached_exposure_rate_air;
    CellCapacity m_capacity;

    explicit Cell(std::vector<observer_ptr<Person<FP>>> persons = {})
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

    ScalarType compute_space_per_person_relative()
    {
    /*
        For every cell in a location we have a transmission factor that is nomalized to m_capacity.volume / m_capacity.persons of
        the location "Home", which is 66. We multiply this rate with the individual size of each cell to obtain a "space per person" factor.
    */
        if (m_capacity.volume != 0) {
            return 66.0 / m_capacity.volume;
        }
        else {
            return 1.0;
        }
    }

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
template<typename FP=double>
class Location
{
public:
    /**
     * @brief Construct a Location of a certain type.
     * @param type The #LocationType.
     * @param index The index of the Location.
     * @param num_cells [Default: 1] The number of Cell%s in which the Location is divided.
     */
    Location(LocationId loc_id, uint32_t num_cells = 1)
        : m_id(loc_id)
        , m_capacity_adapted_transmission_risk(false)
        , m_cells(num_cells)
        , m_required_mask(MaskType::Community)
        , m_npi_active(false)
    {
        assert(num_cells > 0 && "Number of cells has to be larger than 0.");
    }

    Location(LocationType loc_type, uint32_t loc_index, uint32_t num_cells = 1)
        : Location(LocationId{loc_index, loc_type}, num_cells)
    {
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
     * @return Amount of average Infection%s with the virus from the AgeGroup of the transmitter per day.
    */
    ScalarType transmission_contacts_per_day(uint32_t cell_index, VirusVariant virus,
                                             AgeGroup age_receiver) const
    {
        ScalarType prob = 0;
        for (uint32_t age_transmitter = 0; age_transmitter != static_cast<uint32_t>(AgeGroup::Count); ++age_transmitter) {
            prob += m_cells[cell_index].m_cached_exposure_rate_contacts[{virus, static_cast<AgeGroup>(age_transmitter)}] *
                    m_parameters.get<ContactRates>()[{age_receiver, static_cast<AgeGroup>(age_transmitter)}];
        }
        return prob;
    }

    /**
     * @brief Compute the transmission factor for a aerosol transmission of the virus in a Cell.
     * @param[in] cell_index Cell index of the Cell.
     * @param[in] virus VirusVariant of interest.
     * @return Amount of average Infection%s with the virus per day.
    */
    ScalarType transmission_air_per_day(uint32_t cell_index, VirusVariant virus) const
    {
        return m_cells[cell_index].m_cached_exposure_rate_air[{virus}] *
               m_parameters.get<AerosolTransmissionRates>()[{virus}];
    }


    /** 
     * @brief A Person interacts with the population at this Location and may become infected.
     * @param[in, out] person The Person that interacts with the population.
     * @param[in] dt Length of the current Simulation time step.
     * @param[in] global_params Global infection parameters.
     */
    void interact(Person<FP>& person, TimePoint t, TimeSpan dt,
                  const GlobalInfectionParameters<FP>& global_params) const
    {
        // TODO: we need to define what a cell is used for, as the loop may lead to incorrect results for multiple cells
        auto age_receiver          = person.get_age();
        ScalarType mask_protection = person.get_mask_protective_factor(global_params);
        assert(person.get_cells().size() && "Person is in multiple cells. Interact logic is incorrect at the moment.");
        for (auto cell_index :
             person.get_cells()) { // TODO: the logic here is incorrect in case a person is in multiple cells
            std::pair<VirusVariant, ScalarType> local_indiv_trans_prob[static_cast<uint32_t>(VirusVariant::Count)];
            for (uint32_t v = 0; v != static_cast<uint32_t>(VirusVariant::Count); ++v) {
                VirusVariant virus = static_cast<VirusVariant>(v);
                ScalarType local_indiv_trans_prob_v =
                    (std::min(m_parameters.get<MaximumContacts>(),
                              transmission_contacts_per_day(cell_index, virus, age_receiver)) +
                     transmission_air_per_day(cell_index, virus)) *
                    (1 - mask_protection) * dt.days() / days(1).days() * person.get_protection_factor(virus, t);

                local_indiv_trans_prob[v] = std::make_pair(virus, local_indiv_trans_prob_v);
            }
            VirusVariant virus =
                random_transition(VirusVariant::Count, dt,
                                  local_indiv_trans_prob); // use VirusVariant::Count for no virus submission
            if (virus != VirusVariant::Count) {
                person.add_new_infection(
                    mio::abm::Infection<FP>(virus, age_receiver, global_params, t + dt / 2)); // Starting time in first approximation
            }
        }
    }


    /** 
     * @brief Add a Person to the population at this Location.
     * @param[in] person The Person arriving.
     * @param[in] cell_idx [Default: 0] Index of the Cell the Person shall go to.
    */
    void add_person(Person<FP>& p, std::vector<uint32_t> cells = {0})
    {
        m_persons.push_back(&p);
        for (uint32_t cell_idx : cells)
            m_cells[cell_idx].m_persons.push_back(&p);
    }


    /** 
     * @brief Remove a Person from the population of this Location.
     * @param[in] person The Person leaving.
     */
    void remove_person(Person<FP>& p)
    {
        m_persons.erase(std::remove(m_persons.begin(), m_persons.end(), &p), m_persons.end());
        for (auto&& cell : m_cells) {
            cell.m_persons.erase(std::remove(cell.m_persons.begin(), cell.m_persons.end(), &p), cell.m_persons.end());
        }
    }

    /** 
     * @brief Prepare the Location for the next Simulation step.
     * @param[in] t Current TimePoint of the Simulation.
     * @param[in] dt The duration of the Simulation step.
     */
    void cache_exposure_rates(TimePoint t, TimeSpan dt)
    {
        //cache for next step so it stays constant during the step while subpopulations change
        //otherwise we would have to cache all state changes during a step which uses more memory
        for (auto& cell : m_cells) {
            cell.m_cached_exposure_rate_contacts = {{VirusVariant::Count, AgeGroup::Count}, 0.};
            cell.m_cached_exposure_rate_air      = {{VirusVariant::Count}, 0.};
            for (auto&& p : cell.m_persons) {
                if (p->is_infected(t)) {
                    auto inf   = p->get_infection();
                    auto virus = inf.get_virus_variant();
                    auto age   = p->get_age();
                    /* average infectivity over the time step
                 *  to second order accuracy using midpoint rule
                */
                    cell.m_cached_exposure_rate_contacts[{virus, age}] += inf.get_infectivity(t + dt / 2);
                    cell.m_cached_exposure_rate_air[{virus}] += inf.get_infectivity(t + dt / 2);
                }
            }
            if (m_capacity_adapted_transmission_risk) {
                cell.m_cached_exposure_rate_air.array() *= cell.compute_space_per_person_relative();
            }
        }
    }



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
    const std::vector<Cell<FP>>& get_cells() const
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
     * @brief Get the contact exposure rate in the Cell.
     * @param[in] cell_idx Cell index of interest.
     * @return Air exposure rate in the Cell.
     */
    CustomIndexArray<ScalarType, VirusVariant, AgeGroup> get_cached_exposure_rate_contacts(uint32_t cell_idx)
    {
        return m_cells[cell_idx].m_cached_exposure_rate_contacts;
    }

    /**
     * @brief Get the air exposure rate in the Cell.
     * @param[in] cell_idx Cell index of interest.
     * @return Contact exposure rate in the cell.
     */
    CustomIndexArray<ScalarType, VirusVariant> get_cached_exposure_rate_air(uint32_t cell_idx)
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
    CellCapacity get_capacity(uint32_t cell_idx = 0)
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
     * @brief Get the total number of Person%s at the Location.
     * @return Number of Person%s.
     */
    size_t get_number_persons()
    {
        return m_persons.size();
    }


    /**
     * @brief Get the number of Person%s of a particular #InfectionState for all Cell%s.
     * @param[in] t TimePoint of querry.
     * @param[in] state #InfectionState of interest.
     * @return Amount of Person%s of the #InfectionState in all Cell%s.
     */
    size_t get_subpopulation(TimePoint t, InfectionState state) const
    {
        return count_if(m_persons.begin(), m_persons.end(), [&](observer_ptr<Person<FP>> p) {
            return p->get_infection_state(t) == state;
        });
    }


    /**
     * Add a TimePoint to the subpopulations TimeSeries.
     * @param[in] t The TimePoint to be added.
     */
    void store_subpopulations(const TimePoint t)
    {
        m_subpopulations.add_time_point(t.days());
        Eigen::VectorXd subpopulations(Eigen::VectorXd::Zero((size_t)InfectionState::Count));
        for (auto p : m_persons)
            ++subpopulations[(size_t)p->get_infection_state(t)];
        m_subpopulations.get_last_value() = subpopulations;
    }

    /**
     * @brief Initialize the history of subpopulations.
     * @param[in] t The TimePoint of initialization.
     */
    void initialize_subpopulations(TimePoint t)
    {
        if (m_subpopulations.get_num_time_points() == 0) {
            store_subpopulations(t);
        }
        else {
            if (m_subpopulations.get_last_time() != t.days()) { // if not already saved
                store_subpopulations(t);
            }
        }
    }

    /**
     * @brief Get the complete history of subpopulations.
     * @return The TimeSeries of the #InfectionState%s for each TimePoint at the Location.
     */
    const TimeSeries<ScalarType>& get_subpopulations() const
    {
        return m_subpopulations;
    }

private:
    LocationId m_id; ///< Id of the Location including type and index.
    bool m_capacity_adapted_transmission_risk; /**< If true considers the LocationCapacity for the computation of the 
    transmission risk.*/
    LocalInfectionParameters m_parameters; ///< Infection parameters for the Location.
    std::vector<observer_ptr<Person<FP>>> m_persons{}; ///< A vector of all Person%s at the Location.
    TimeSeries<ScalarType> m_subpopulations{Eigen::Index(
        InfectionState::Count)}; ///< A TimeSeries of the #InfectionState%s for each TimePoint at the Location.
    std::vector<Cell<FP>> m_cells{}; ///< A vector of all Cell%s that the Location is divided in.
    MaskType m_required_mask; ///< Least secure type of Mask that is needed to enter the Location.
    bool m_npi_active; ///< If true requires e.g. Mask%s to enter the Location.
};

} // namespace abm
} // namespace mio

#endif
