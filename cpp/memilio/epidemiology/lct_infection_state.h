/*
* Copyright (C) 2020-2025 MEmilio
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

#include "memilio/config.h"
#include "memilio/math/eigen.h"

#include <array>

namespace mio
{
/**
 * @brief Provides the functionality to be able to work with subcompartments in an LCT model.

 * This class just stores the number of subcompartments for each InfectionState and not the number of individuals in
 * each subcompartment.
 *
 * @tparam InfectionStates An enum class that defines the basic infection states.
 * @tparam Ns Number of subcompartments for each infection state defined in InfectionState.
 *      The number of given template arguments must be equal to the entry Count from InfectionStates.
 */
template <typename FP, class InfectionStates, size_t... Ns>
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

    /**
     * @brief Cumulates a vector with the number of individuals in each subcompartment (with subcompartments
     *  according to the LctInfectionState) to produce a Vector that divides the population only into the infection
     *  states defined in InfectionStates.
     *
     * @param[in] subcompartments Vector with number of individuals in each subcompartment.
     *  The size of the vector has to match the LctInfectionState.
     * @return Vector with accumulated values for the InfectionStates.
     */
    static Eigen::VectorX<FP> calculate_compartments(const Eigen::VectorX<FP>& subcompartments)
    {
        assert(subcompartments.rows() == Count);

        Eigen::VectorX<FP> compartments((Eigen::Index)InfectionState::Count);
        // Use segment of the vector subcompartments of each InfectionState and sum up the values of subcompartments.
        for (int i = 0; i < (Eigen::Index)InfectionState::Count; i++) {
            InfectionState State = static_cast<InfectionState>(i);
            // first index of first subcompartment:
            size_t index = 0;
            for (size_t j = 0; j < (size_t)(State); j++) {
                index = index + m_subcompartment_numbers[j];
            }
            // number of subcompartments:
            size_t num_subcomp = m_subcompartment_numbers[(size_t)State];
            compartments[i]    = subcompartments.segment(index, num_subcomp).sum();
        }

        return compartments;
    }

    static constexpr size_t Count{(... + Ns)};

private:
    static constexpr const std::array<size_t, sizeof...(Ns)> m_subcompartment_numbers{
        Ns...}; ///< Vector which defines the number of subcompartments for each infection state of InfectionState.
};

} // namespace mio

#endif
