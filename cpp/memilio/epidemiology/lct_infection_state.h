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
#ifndef MIO_EPI_LCT_INFECTION_STATE_H
#define MIO_EPI_LCT_INFECTION_STATE_H

#include <array>

namespace mio
{
/**
 * @brief Provides the functionality to be able to work with subcompartments in an LCT model.
 *
 * @tparam InfectionStates An enum class that defines the basic infection states.
 * @tparam Ns Number of subcompartments for each infection state defined in InfectionState. 
 *      The number of given template arguments must be equal to the entry Count from InfectionState.
 */
template <class InfectionStates, size_t... Ns>
class LctInfectionState
{
public:
    using InfectionState = InfectionStates;
    static_assert((size_t)InfectionState::Count == sizeof...(Ns),
                  "The number of the size_t's provided as template parameters must be "
                  "the same as the entry Count of InfectionState.");

    static_assert(((Ns > 0) && ...), "The number of subcompartments must be at least 1.");

    /**
     * @brief Gets the number of subcompartments in an infection state.
     *
     * @tparam State Infection state for which the number of subcompartments should be returned.   
     * @return Number of subcompartments for State. Returned value is always at least one.
     */
    template <InfectionState State>
    static constexpr size_t get_num_subcompartments()
    {
        static_assert(State < InfectionState::Count, "State must be a a valid InfectionState.");
        return m_subcompartment_numbers[(size_t)State];
    }

    /**
     * @brief Gets the index of the first subcompartment of an infection state.
     *
     * In a simulation, the number of individuals in the subcompartments are stored in vectors. 
     * Accordingly, the index of the first subcompartment of State in such a vector is returned.
     * @tparam State: Infection state for which the index should be returned.    
     * @return Index of the first subcompartment for a vector with one entry per subcompartment. 
     *      Returned value is always non-negative.
     */
    template <InfectionState State>
    static constexpr size_t get_first_index()
    {
        static_assert(State < InfectionState::Count, "State must be a a valid InfectionState.");
        size_t index = 0;
        for (size_t i = 0; i < (size_t)(State); i++) {
            index = index + m_subcompartment_numbers[i];
        }
        return index;
    }

    static constexpr size_t Count{(... + Ns)};

private:
    static constexpr const std::array<size_t, sizeof...(Ns)> m_subcompartment_numbers{
        Ns...}; ///< Vector which defines the number of subcompartments for each infection state of InfectionState.
};

} // namespace mio

#endif