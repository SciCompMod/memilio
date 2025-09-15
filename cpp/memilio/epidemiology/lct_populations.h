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
#ifndef MIO_EPI_LCT_POPULATIONS_H
#define MIO_EPI_LCT_POPULATIONS_H

#include "boost/type_traits/make_void.hpp"
#include "memilio/config.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/math/eigen.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/utils/type_list.h"
#include "memilio/utils/metaprogramming.h"

namespace mio
{

/**
 * @brief A class template for compartment populations of LCT models.
 *
 * Populations can be split up into different categories, e.g. by age group, yearly income group, gender etc.
 * In LCT models, we want to use different numbers of subcompartments, i.e. different LctStates,
 * for each group of a category.
 * (Therefore, we can't use the normal Populations class because it expects the same InfectionStates for each group.)
 *
 * This template must be instantiated with one LctState for each group of a category.
 * The purpose of the LctStates is to define the number of subcompartments for each InfectionState.
 * The number of LctStates also determines the number of groups.
 * If you want to use more than one category, e.g. age and gender, you have to define num_age_groups * num_genders
 * LctStates, because the number of subcompartments can be different
 * even for (female, A05-A14) and (female, A80+) or (male, A05-A14).
 *
 * The class created from this template contains a "flat array" of compartment populations and some functions for
 * retrieving or setting the populations. The order in the flat array is: First, all (sub-)compartments of the
 * first group, afterwards all (sub-)compartments of the second group and so on.
 *
 */

template <typename FP, class... LctStates>
class LctPopulations
{
public:
    using Type                         = UncertainValue<FP>;
    using InternalArrayType            = Eigen::Array<Type, Eigen::Dynamic, 1>;
    using LctStatesGroups              = TypeList<LctStates...>;
    static size_t constexpr num_groups = sizeof...(LctStates); ///< Number of groups.
    static_assert(num_groups >= 1, "The number of LctStates provided should be at least one.");

    /// @brief Default constructor.
    LctPopulations()
    {
        set_count();
        m_y = InternalArrayType::Constant(m_count, UncertainValue<FP>(0.0));
    }

    /**
     * @brief get_num_compartments Returns the number of compartments.
     * @return Number of compartments which equals the flat array size.
     */
    size_t get_num_compartments() const
    {
        return m_count;
    }

    /**
     * @brief Returns a reference to the internally stored flat array.
     * @return Const reference to the InternalArrayType instance.
     */
    auto const& array() const
    {
        return m_y;
    }
    auto& array()
    {
        return m_y;
    }

    /**
     * @brief Returns the entry of the array given a flat index.
     * @param index A flat index.
     * @return The value of the internal array at the index.
     */
    Type& operator[](size_t index)
    {
        assert(index < m_count);
        return m_y[index];
    }

    /**
    * @brief Gets the first index of a group in the flat array.
    * @tparam group The group for which the index should be returned.
    * @return The index of the first entry of group in the flat array.
    */
    template <size_t Group>
    size_t get_first_index_of_group() const
    {
        static_assert((Group < num_groups) && (Group >= 0), "The template parameter Group should be valid.");
        if constexpr (Group == 0) {
            return 0;
        }
        else {
            return get_first_index_of_group<Group - 1>() + type_at_index_t<Group - 1, LctStatesGroups>::Count;
        }
    }
    /**
     * @brief Returns an Eigen copy of the vector of populations.
     * This can be used as initial conditions for the ODE solver.
     * @return Eigen::VectorXd of populations.
     */
    inline Eigen::VectorX<FP> get_compartments() const
    {
        return m_y.array().template cast<FP>();
    }

    /**
     * @brief Returns the total population of a group.
     * @tparam group The group for which the total population should be calculated.
     * @return Total population of the group.
     */
    template <size_t Group>
    FP get_group_total() const
    {
        return m_y.array()
            .template cast<FP>()
            .segment(get_first_index_of_group<Group>(), type_at_index_t<Group, LctStatesGroups>::Count)
            .sum();
    }

    /**
     * @brief Returns the total population of all compartments and groups.
     * @return Total population.
     */
    FP get_total() const
    {
        return m_y.array().template cast<FP>().sum();
    }

    /**
     * @brief Checks whether all compartments have non-negative values.
     * This function can be used to prevent slightly negative function values in compartment sizes that are produced
     * due to rounding errors if, e.g., population sizes were computed in a complex way.
     *
     * Attention: This function should be used with care. It can not and will not set model parameters and
     *            compartments to meaningful values. In most cases it is preferable to use check_constraints,
     *            and correct values manually before proceeding with the simulation.
     *            The main usage for apply_constraints is in automated tests using random values for initialization.
     *
     * @return Returns true if one (or more) constraint(s) were corrected, otherwise false.
    */
    bool apply_constraints()
    {
        bool corrected = false;
        for (int i = 0; i < m_y.array().size(); i++) {
            if (m_y.array()[i] < 0.0) {
                log_warning("Constraint check: Compartment size {:d} changed from {:.4f} to {:d}", i, m_y.array()[i],
                            0);
                m_y.array()[i] = 0.;
                corrected      = true;
            }
        }
        return corrected;
    }

    /**
     * @brief Checks whether all compartments have non-negative values and logs an error if constraint is not satisfied.
     * @return Returns true if one or more constraints are not satisfied, false otherwise.
     */
    bool check_constraints() const
    {
        if ((m_y.array() < 0.0).any()) {
            log_error("Constraint check: At least one compartment size is smaller {}.", 0);
            return true;
        }
        return false;
    }

private:
    /**
     * @brief Sets recursively the total number of (sub-)compartments over all groups.
     * The number also corresponds to the size of the internal vector.
     */
    template <size_t Group = 0>
    void set_count()
    {
        if constexpr (Group == 0) {
            m_count = 0;
        }
        if constexpr (Group < num_groups) {
            m_count += type_at_index_t<Group, LctStatesGroups>::Count;
            set_count<Group + 1>();
        }
    }

    size_t m_count; //< Number of groups stored.
    InternalArrayType m_y{}; //< An array containing the number of people in the groups.
};

} // namespace mio

#endif // MIO_EPI_LCT_POPULATIONS_H
