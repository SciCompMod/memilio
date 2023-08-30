/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#ifndef COMPARTMENTALMODEL_H
#define COMPARTMENTALMODEL_H

#include "memilio/config.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/metaprogramming.h"
#include <vector>
#include <functional>

#define USE_DERIV_FUNC 1

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
 * @brief check whether a check_constraints function exists
 * @tparam The type to check for the existence of the member function
 */
template <class T>
using has_check_constraints_member_function = is_expression_valid<details::check_constraints_expr_t, T>;

/**
 * @brief check whether a apply_constraints function exists
 * @tparam The type to check for the existence of the member function
 */
template <class T>
using has_apply_constraints_member_function = is_expression_valid<details::apply_constraints_expr_t, T>;

/**
 * @brief CompartmentalModel is a template for a compartmental model for an
 * array of initial populations and a parameter set
 *
 * The Populations must be a concrete class derived from the Populations template,
 * i.e. a multi-dimensional array of compartment populations where each dimension
 * corresponds to a category
 *
 * The ParameterSet must be a concrete class derived form the ParameterSet template,
 * i.e. a compile-time map of parameters used by the model. These can be referenced
 * when defining the flows between compartments and they can be used for parameter
 * studies
 *
 */
template <class Comp, class Pop, class Params>
struct CompartmentalModel {
public:
    using Compartments = Comp;
    using Populations  = Pop;
    using ParameterSet = Params;

    // The flow function takes a set of parameters, the current time t and the
    // snapshot y of all population sizes at time t, represented as a flat array and returns a scalar value
    // that represents a flow going from one compartment to another.
    using FlowFunction = std::function<ScalarType(ParameterSet const& p, Eigen::Ref<const Eigen::VectorXd> pop,
                                                  Eigen::Ref<const Eigen::VectorXd> y, double t)>;

    // A flow is a tuple of a from-index corresponding to the departing compartment, a to-index
    // corresponding to the receiving compartment and a FlowFunction. The value returned by the flow
    // function will be subtracted from the time derivative of the populations at the flat index corresponding
    // to the from-compartment, and added to the time derivative of the populations at the flat index
    // corresponding to the to-compartment.
    using Flow = std::tuple<typename Populations::Index, typename Populations::Index, FlowFunction>;

    /**
     * @brief CompartmentalModel default constructor
     */
    CompartmentalModel(Populations const& po, ParameterSet const& pa)
        : populations{std::move(po)}
        , parameters{pa}
    {
    }

    CompartmentalModel(const CompartmentalModel&)            = default;
    CompartmentalModel(CompartmentalModel&&)                 = default;
    CompartmentalModel& operator=(const CompartmentalModel&) = default;
    CompartmentalModel& operator=(CompartmentalModel&&)      = default;
    virtual ~CompartmentalModel()                            = default;

    /**
     * @brief add_flow defines a flow from compartment A to another compartment B
     * @param from is the index of the departure compartment A
     * @param to is the index of the receiving compartment B
     * @param f is a function defining the flow given a set of parameters, the current time t and the
     * snapshot y of all population sizes at time t, represented as a flat array
     */
    void add_flow(typename Populations::Index from, typename Populations::Index to, FlowFunction f)
    {
        flows.push_back(Flow(from, to, f));
    }

#if USE_DERIV_FUNC
    //REMARK: Not pure virtual for easier java/python bindings
    virtual void get_derivatives(Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<const Eigen::VectorXd> /*y*/,
                                 double /*t*/, Eigen::Ref<Eigen::VectorXd> /*dydt*/) const {};
#endif // USE_DERIV_FUNC

    /**
     * @brief eval_right_hand_side evaulates the right-hand-side f of the ODE dydt = f(y, t)
     *
     * The heart of the compartmental model is a first order ODE dydt = f(y,t), where y is a flat
     * representation of all the compartmental populations at time t. This function evaluates the
     * right-hand-side f of the ODE from the intercompartmental flows. It can be used in an ODE
     * solver
     *
     * The distinction between pop and y is only for the case of mobility.
     * If we have mobility, we want to evaluate the evolution of infection states for a small group of travellers (y)
     * while they are in any population (pop). It is important that pop > y always applies.
     *
     * If we consider a simulation without mobility, the function is called with
     * model.eval_right_hand_side(y, y, t, dydt)
     *
     * @param pop the current state of the population in the geographic unit we are considering
     * @param y the current state of the model (or a subpopulation) as a flat array
     * @param t the current time
     * @param dydt a reference to the calculated output
     */
    void eval_right_hand_side(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                              Eigen::Ref<Eigen::VectorXd> dydt) const
    {
        dydt.setZero();

#if USE_DERIV_FUNC
        this->get_derivatives(pop, y, t, dydt);
#else // USE_DERIV_FUNC
        for (auto& flow : flows) {
            ScalarType f = std::get<2>(flow)(parameters, pop, y, t);
            dydt[call(Populations::get_flat_index, std::get<0>(flow))] -= f;
            dydt[call(Populations::get_flat_index, std::get<1>(flow))] += f;
        }
#endif // USE_DERIV_FUNC
    }

    /**
     * @brief get_initial_values returns the initial values for the compartmental populations.
     * This can be used as initial conditions in an ODE solver
     * @return the initial populatoins
     */
    Eigen::VectorXd get_initial_values() const
    {
        return populations.get_compartments();
    }

    // TODO: if constexpr as soon as we open for C++17
    template <typename T = ParameterSet>
    std::enable_if_t<has_apply_constraints_member_function<T>::value> apply_constraints()
    {
        populations.apply_constraints();
        parameters.apply_constraints();
    }

    template <typename T = ParameterSet>
    std::enable_if_t<!has_apply_constraints_member_function<T>::value> apply_constraints()
    {
        populations.apply_constraints();
    }

    template <typename T = ParameterSet>
    std::enable_if_t<has_check_constraints_member_function<T>::value> check_constraints() const
    {
        populations.check_constraints();
        parameters.check_constraints();
    }

    template <typename T = ParameterSet>
    std::enable_if_t<!has_check_constraints_member_function<T>::value> check_constraints() const
    {
        populations.check_constraints();
    }

    Populations populations{};
    ParameterSet parameters{};

private:
    std::vector<Flow> flows{};
};

/**
 * detect the eval_right_hand_side member function of a compartment model.
 * If the eval_right_hand_side member function exists in the type M, this template when instatiated
 * will be equal to the return type of the function.
 * Otherwise the template is invalid.
 * @tparam M a type that has a eval_right_hand_side member function, e.g. a compartment model type.
 */
template <class M>
using eval_right_hand_side_expr_t = decltype(std::declval<const M&>().eval_right_hand_side(
    std::declval<Eigen::Ref<const Eigen::VectorXd>>(), std::declval<Eigen::Ref<const Eigen::VectorXd>>(),
    std::declval<double>(), std::declval<Eigen::Ref<Eigen::VectorXd>>()));

/**
 * detect the get_initial_values member function of a compartment model.
 * If the detect_initial_values member function exists in the type M, this template when instatiated
 * will be equal to the return type of the function.
 * Otherwise the template is invalid.
 * @tparam M a type that has a get_initial_values member function, e.g. a compartment model type.
 */
template <class M>
using get_initial_values_expr_t =
    decltype(std::declval<Eigen::VectorXd&>() = std::declval<const M&>().get_initial_values());

/**
 * Template meta function to check if a type is a valid compartment model. 
 * Defines a static constant of name `value`. 
 * The constant `value` will be equal to true if M is a valid compartment model type.
 * Otherwise, `value` will be equal to false.
 * @tparam Sim a type that may or may not be a compartment model.
 */
template <class M>
using is_compartment_model = std::integral_constant<bool, (is_expression_valid<eval_right_hand_side_expr_t, M>::value &&
                                                           is_expression_valid<get_initial_values_expr_t, M>::value)>;

} // namespace mio

#endif // COMPARTMENTALMODEL_H
