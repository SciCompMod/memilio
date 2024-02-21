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
#include "memilio/math/eigen.h"
#include <array>

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
 * @brief Provides the functionality to be able to work with subcompartments in an LCT model.
 *
 * @tparam InfectionStateBase An enum class that defines the basic infection states.
 * @tparam Ns Number of subcompartments for each infection state defined in InfectionStateBase. 
 *      The number of given template arguments must be equal to the entry Count from InfectionStateBase.
 */
template <class InfectionStateBase, unsigned int... Ns>
class InfectionState
{
    using Base = InfectionStateBase;
    static_assert((unsigned int)Base::Count == sizeof...(Ns),
                  "The number of integers provided as template parameters must be "
                  "the same as the entry Count of InfectionStateBase.");

    static_assert(((Ns > 0) && ...), "The number of subcompartments must be at least 1.");

public:
    /**
     * @brief Gets the number of subcompartments in an infection state.
     *
     * @param[in] infectionstatebase Infection state for which the number of subcompartments should be returned.   
     * @return Number of subcompartments for infectionstatebase.
     */
    template <Base infectionstatebase>
    static constexpr unsigned int get_number()
    {
        return m_subcompartment_numbers[(int)infectionstatebase];
    }

    /**
     * @brief Gets the index of the first subcompartment in an vector with all subcompartments for template an infection state.
     *
     * In a simulation, the number of individuals in the subcompartments are stored in vectors. 
     * Accordingly, the index in such a vector of the first subcompartment of an infection state is given.
     * @param[in] infectionstatebase Infection state for which the index should be returned.    
     * @return Index of the first subcompartment for a vector with one entry per subcompartment.
     */
    template <Base infectionstatebase>
    static constexpr unsigned int get_firstindex()
    {
        unsigned int index = 0;
        for (int i = 0; i < (int)(infectionstatebase); i++) {
            index = index + m_subcompartment_numbers[i];
        }
        return index;
    }

    static constexpr unsigned int Count{(... + Ns)};

private:
    static constexpr const std::array<unsigned int, sizeof...(Ns)> m_subcompartment_numbers{
        Ns...}; ///< Vector which defines the number of subcompartments for each infection state of InfectionStateBase.
};

} // namespace lsecir
} // namespace mio

#endif