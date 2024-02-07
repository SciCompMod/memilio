/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Lena Ploetzke
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

#ifndef LCTSECIR_INFECTIONSTATE_H
#define LCTSECIR_INFECTIONSTATE_H

#include "memilio/utils/logging.h"
#include <vector>

namespace mio
{
namespace lsecir
{

/**
 * @brief The InfectionStateBase enum describes the possible basic
 * categories for the infection state of persons
 */
enum class InfectionStateBase
{
    Susceptible        = 0,
    Exposed            = 1,
    InfectedNoSymptoms = 2,
    InfectedSymptoms   = 3,
    InfectedSevere     = 4,
    InfectedCritical   = 5,
    Recovered          = 6,
    Dead               = 7,
    Count              = 8
};

/**
 * @brief The InfectionTransition enum describes the possible
 * transitions of the infectious state of persons.
 */
enum class InfectionTransition
{
    SusceptibleToExposed                 = 0,
    ExposedToInfectedNoSymptoms          = 1,
    InfectedNoSymptomsToInfectedSymptoms = 2,
    InfectedNoSymptomsToRecovered        = 3,
    InfectedSymptomsToInfectedSevere     = 4,
    InfectedSymptomsToRecovered          = 5,
    InfectedSevereToInfectedCritical     = 6,
    InfectedSevereToRecovered            = 7,
    InfectedCriticalToDead               = 8,
    InfectedCriticalToRecovered          = 9,
    Count                                = 10
};

class InfectionState
{
public:
    /**
     * @brief Constructor for the InfectionState class.
     *
     * InfectionState class defines the possible InfectionState%s with the number of Subcompartments for the LCT model.
     * With the default constructor, the class is defined without subcompartments, i.e. only the subdivision in InfectionStateBase 
     * is used.
     */
    InfectionState()
        : m_subcompartment_numbers((int)InfectionStateBase::Count, 1)
        , m_subcompartment_indexfirst((int)InfectionStateBase::Count, 1)
    {
        set_compartment_index();
    }

    /**
     * @brief Constructor for the InfectionState class.
     *
     * InfectionState class defines the possible InfectionState%s with the number of Subcompartments for the LCT model.
     * @param[in] subcompartment_numbers Vector which defines the number of Subcompartments for each infection state of InfectionStateBase.       
     */
    InfectionState(std::vector<int> subcompartment_numbers)
        : m_subcompartment_numbers(std::move(subcompartment_numbers))
        , m_subcompartment_indexfirst((int)InfectionStateBase::Count, 1)
    {
        check_constraints();
        set_compartment_index();
    }

    /**
     * @brief Setter for the number of Subcompartments.
     *
     * The number of Subcompartments is only updated if the vector is valid.
     * @param[in] subcompartment_numbers Vector which defines the number of Subcompartments for each infection state of InfectionStateBase. 
     * @return Returns true if the vector is not valid, otherwise false.      
     */
    bool set_subcompartment_numbers(std::vector<int> subcompartment_numbers)
    {
        std::vector<int> copy_m_subcompartment_numbers(m_subcompartment_numbers);
        m_subcompartment_numbers = std::move(subcompartment_numbers);
        if (check_constraints()) {
            // Case where the vector is not valid.
            m_subcompartment_numbers = copy_m_subcompartment_numbers;
            return true;
        }
        else {
            set_compartment_index();
            return false;
        }
    }

    /**
     * @brief Gets the number of subcompartments in an infection state.
     *
     * @param[in] infectionstatebase Infection state for which the number of subcompartments should be returned.   
     * @return Number of Subcompartments for infectionstatebase.
     */
    int get_number(InfectionStateBase infectionstatebase) const
    {
        return m_subcompartment_numbers[(int)infectionstatebase];
    }

    /**
     * @brief Gets the number of subcompartments in an infection state.
     *
     * @param[in] infectionstatebase Index of an infection state for which the number of subcompartments should be returned.   
     * If the index does not match an infectionstate, the return value will be -1.
     * @return Number of Subcompartments for infectionstatebase or -1.
     */
    int get_number(int infectionstatebaseindex) const
    {
        if ((0 <= infectionstatebaseindex) && (infectionstatebaseindex < (int)InfectionStateBase::Count)) {
            return m_subcompartment_numbers[infectionstatebaseindex];
        }
        else {
            return -1;
        }
    }

    /**
     * @brief Gets the index of the first subcompartment in an vector with all subcompartments for an infection state.
     *
     * In a simulation, the number of individuals in the subcompartments are stored in vectors. 
     * Accordingly, the index in such a vector of the first subcompartment of an infection state is given.
     * @param[in] infectionstatebase Infection state for which the index should be returned.    
     * @return Index of the first Subcompartment for a vector with one entry per subcompartment.
     */
    int get_firstindex(InfectionStateBase infectionstatebase) const
    {
        return m_subcompartment_indexfirst[(int)infectionstatebase];
    }

    /**
     * @brief Gets the index of the first subcompartment in an vector with all subcompartments for an infection state.
     *
     * In a simulation, the number of individuals in the subcompartments are stored in vectors. 
     * Accordingly, the index in such a vector of the first subcompartment of an infection state is given.
     * @param[in] infectionstatebase Index of an infection state for which the index of a vector should be returned.   
     * If the index does not match an infectionstate, the return value will be -1.  
     * @return Index of the first Subcompartment for a vector with one entry per subcompartment or -1.
     */
    int get_firstindex(int infectionstatebaseindex) const
    {
        if ((0 <= infectionstatebaseindex) && (infectionstatebaseindex < (int)InfectionStateBase::Count)) {
            return m_subcompartment_indexfirst[infectionstatebaseindex];
        }
        else {
            return -1;
        }
    }

    /**
     * @brief Gets the total number of (sub-)compartments of infection states.
     */
    int get_count() const
    {
        return m_count;
    }

    /**
     * @brief Checks constraints on infection states.
     *
     * @return Returns true if one (or more) constraint(s) are not satisfied, otherwise false.
     */
    bool check_constraints() const
    {
        if (!(m_subcompartment_numbers.size() == (int)InfectionStateBase::Count)) {
            log_error("Vector for number of subcompartments has the wrong size.");
            return true;
        }
        if (!(m_subcompartment_numbers[(int)InfectionStateBase::Susceptible] == 1)) {
            log_error("Susceptible compartment can not have subcompartments.");
            return true;
        }
        if (!(m_subcompartment_numbers[(int)InfectionStateBase::Recovered] == 1)) {
            log_error("Recovered compartment can not have subcompartments.");
            return true;
        }
        if (!(m_subcompartment_numbers[(int)InfectionStateBase::Dead] == 1)) {
            log_error("Dead compartment can not have Subcompartments.");
            return true;
        }
        for (int i = 0; i < (int)InfectionStateBase::Count; ++i) {
            if (m_subcompartment_numbers[i] < 1) {
                log_error("All compartments should have at least one Subcompartment.");
                return true;
            }
        }
        return false;
    }

private:
    /**
     * @brief Calculates Index of the first Subcompartment for a vector with one entry per subcompartment.
     *
     * Therefore the vector with number of subcompartments per infection state is used.
     */
    void set_compartment_index()
    {
        int index = 0;
        for (int i = 0; i < (int)(InfectionStateBase::Count); i++) {
            m_subcompartment_indexfirst[i] = index;
            index                          = index + m_subcompartment_numbers[i];
        }
        m_count = index;
    }

    std::vector<int>
        m_subcompartment_numbers; ///< Vector which defines the number of Subcompartments for each infection state of InfectionStateBase.
    std::vector<int>
        m_subcompartment_indexfirst; ///< Vector with Indexes for all infection states of the first Subcompartment for a vector with one entry per subcompartment.
    int m_count; ///< Total number of (sub-)compartments of infection states.
};

} // namespace lsecir
} // namespace mio

#endif