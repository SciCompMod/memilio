/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Jan Kleinert, Daniel Abele
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
#ifndef MIO_COMPARTMENTALMODEL_H
#define MIO_COMPARTMENTALMODEL_H

#include "memilio/config.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/metaprogramming.h"
#include <cstddef>
#include <type_traits>
#include <vector>
#include <functional>

namespace mio
{

namespace details
{
//helpers for check_constraint
template <class T>
using check_constraints_expr_t = decltype(std::declval<T>().check_constraints());

//helpers for apply_constraints
template <class T>
using apply_constraints_expr_t = decltype(std::declval<T>().apply_constraints());

} //namespace details

/**
 * @brief Check whether a check_constraints function exists.
 * @tparam The type to check for the existence of the member function.
 */
template <class T>
using has_check_constraints_member_function = is_expression_valid<details::check_constraints_expr_t, T>;

/**
 * @brief Check whether a apply_constraints function exists.
 * @tparam The type to check for the existence of the member function.
 */
template <class T>
using has_apply_constraints_member_function = is_expression_valid<details::apply_constraints_expr_t, T>;

/**
 * @brief CompartmentalModel is a template for a compartmental model for an
 * array of initial populations and a parameter set.
 * @tparam FP floating point type, e.g., double.
 *
 * The Populations must be a concrete class derived from the Populations template,
 * i.e. a multi-dimensional array of compartment populations where each dimension
 * corresponds to a category.
 *
 * The ParameterSet must be a concrete class derived form the ParameterSet template,
 * i.e. a compile-time map of parameters used by the model. These can be referenced
 * when defining the flows between compartments and they can be used for parameter
 * studies.
 *
 */
template <typename FP, class Comp, class Pop, class Params>
struct CompartmentalModel {
public:
    using Compartments = Comp;
    using Populations  = Pop;
    using ParameterSet = Params;
    /**
     * @brief CompartmentalModel default constructor.
     */
    CompartmentalModel(Populations const& po, ParameterSet const& pa)
        : populations{std::move(po)}
        , parameters{pa}
    {
    }

    CompartmentalModel(const CompartmentalModel&) = default;
    CompartmentalModel(CompartmentalModel&&)      = default;
    CompartmentalModel& operator=(const CompartmentalModel&) = default;
    CompartmentalModel& operator=(CompartmentalModel&&) = default;
    virtual ~CompartmentalModel()                       = default;

    // REMARK: Not pure virtual for easier java/python bindings.
    virtual void get_derivatives(Eigen::Ref<const Eigen::VectorX<FP>> /*pop*/,
                                 Eigen::Ref<const Eigen::VectorX<FP>> /*y*/, FP /*t*/,
                                 Eigen::Ref<Eigen::VectorX<FP>> /*dydt*/) const {};

    /**
     * @brief This function evaluates the right-hand-side f of the ODE dydt = f(y, t).
     *
     * The heart of the compartmental model is a first order ODE dydt = f(y,t), where y is a flat
     * representation of all the compartmental populations at time t. This function evaluates the
     * right-hand-side f of the ODE from the intercompartmental flows. It can be used in an ODE
     * solver.
     *
     * The distinction between pop and y is only for the case of mobility.
     * If we have mobility, we want to evaluate the evolution of infection states for a small group of travellers (y)
     * while they are in any population (pop). It is important that pop > y always applies.
     *
     * If we consider a simulation without mobility, the function is called with
     * model.eval_right_hand_side(y, y, t, dydt).
     *
     * @param[in] pop The current state of the population in the geographic unit we are considering.
     * @param[in] y The current state of the model (or a subpopulation) as a flat array.
     * @param[in] t The current time.
     * @param[out] dydt A reference to the calculated output.
     */
    void eval_right_hand_side(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                              Eigen::Ref<Eigen::VectorX<FP>> dydt) const
    {
        dydt.setZero();
        this->get_derivatives(pop, y, t, dydt);
    }

    /**
     * @brief Get the initial conditions for the ODE dydt = f(y, t).
     * See eval_right_hand_side for more detail.
     * @return Current value of model populations as a flat vector.
     */
    Eigen::VectorX<FP> get_initial_values() const
    {
        return populations.get_compartments();
    }

    /**
     * @brief Checks whether the model satisfies all constraints. If not, it changes values to suffice their constraints.
     *
     * Attention: This function should be used with care. It can not and will not set model parameters and 
     *            compartments to meaningful values. In most cases it is preferable to use check_constraints,
     *            and correct values manually before proceeding with the simulation.
     *            The main usage for apply_constraints is in automated tests using random values for initialization.
     *
     * @return Returns true if one or more constraints were corrected, false otherwise. 
     */
    bool apply_constraints()
    {
        if constexpr (has_apply_constraints_member_function<ParameterSet>::value) {
            // use bitwise instead of logical or, so that both apply_constraint functions are called
            return ((int)parameters.apply_constraints() | (int)populations.apply_constraints());
        }
        else {
            return populations.check_constraints();
        }
    }

    /**
     * @brief Checks that the model satisfies all constraints (e.g. parameter or population constraints).
     * @return Returns true if one or more constraints are not satisfied, false otherwise.
     */
    bool check_constraints() const
    {
        if constexpr (has_check_constraints_member_function<ParameterSet>::value) {
            return (parameters.check_constraints() || populations.check_constraints());
        }
        else {
            return populations.check_constraints();
        }
    }

    Populations populations{};
    ParameterSet parameters{};
};

/**
 * Detect the eval_right_hand_side member function of a compartment model.
 * If the eval_right_hand_side member function exists in the type M, this template when instatiated
 * will be equal to the return type of the function.
 * Otherwise the template is invalid.
 * @tparam FP, floating point type, e.g., double.
 * @tparam M a type that has a eval_right_hand_side member function, e.g. a compartment model type.
 */
template <typename FP, class M>
using eval_right_hand_side_expr_t = decltype(std::declval<const M&>().eval_right_hand_side(
    std::declval<Eigen::Ref<const Eigen::VectorX<FP>>>(), std::declval<Eigen::Ref<const Eigen::VectorX<FP>>>(),
    std::declval<FP>(), std::declval<Eigen::Ref<Eigen::VectorX<FP>>>()));

/**
 * Detect the get_initial_values member function of a compartment model.
 * If the detect_initial_values member function exists in the type M, this template when instatiated
 * will be equal to the return type of the function.
 * Otherwise the template is invalid.
 * @tparam FP, floating point type, e.g., double.
 * @tparam M a type that has a get_initial_values member function, e.g. a compartment model type.
 */
template <typename FP, class M>
using get_initial_values_expr_t =
    decltype(std::declval<Eigen::VectorX<FP>&>() = std::declval<const M&>().get_initial_values());

/**
 * Template meta function to check if a type is a valid compartment model. 
 * Defines a static constant of name `value`. 
 * The constant `value` will be equal to true if M is a valid compartment model type.
 * Otherwise, `value` will be equal to false.
 * @tparam FP, floating point type, e.g., double.
 * @tparam M a type that may or may not be a compartment model.
 */
template <typename FP, class M>
using is_compartment_model =
    std::integral_constant<bool, (is_expression_valid<eval_right_hand_side_expr_t, FP, M>::value &&
                                  is_expression_valid<get_initial_values_expr_t, FP, M>::value)>;

} // namespace mio

#endif // MIO_COMPARTMENTALMODEL_H
