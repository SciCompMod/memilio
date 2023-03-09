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
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/flow_chart.h"
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
template <class Comp, class Pop, class Params, class FlowChart = void>
struct CompartmentalModel {
public:
    using Compartments = Comp;
    using Populations  = Pop;
    using ParameterSet = Params;

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

    //REMARK: Not pure virtual for easier java/python bindings
    virtual void get_derivatives(Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<const Eigen::VectorXd> /*y*/,
                                 double /*t*/, Eigen::Ref<Eigen::VectorXd> /*dydt*/) const {};
    /**
     * @brief eval_right_hand_side evaulates the right-hand-side f of the ODE dydt = f(y, t)
     *
     * The heart of the compartmental model is a first order ODE dydt = f(y,t), where y is a flat
     * representation of all the compartmental populations at time t. This function evaluates the
     * right-hand-side f of the ODE from the intercompartmental flows. It can be used in an ODE
     * solver
     *
     * @param y the current state of the model as a flat array
     * @param t the current time
     * @param dydt a reference to the calculated output
     */
    void eval_right_hand_side(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                              Eigen::Ref<Eigen::VectorXd> dydt) const
    {
        dydt.setZero();
        this->get_derivatives(pop, y, t, dydt);
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
};

/**
 * @brief CompartmentalModel with flows. 
 *
 * Inherits from CompartmentalModel extending it with a FlowChart (the optional template parameter). 
 * @see CompartmentalModel
 *
 * Instead of directly computing the right-hand-side with get_derivatives, first the flows defined by FlowChart are
 * computed, and then added to/subtracted from their source/target compartment.
 */
template <class Comp, class Pop, class Params, class... Flows>
struct CompartmentalModel<Comp, Pop, Params, FlowChart<Flows...>> : public CompartmentalModel<Comp, Pop, Params> {
    using PopIndex = typename Pop::Index;
    // PopIndex without Comp, used as argument type for get_flow_index
    using FlowIndex = filtered_index_t<Comp, Index, PopIndex>;
    // Enforce that Comp is a unique Category of PopIndex, since we use source/target of Flows to provide Comp indices
    static_assert(FlowIndex::size == Pop::Index::size - 1,
                  "Compartments must be used exactly once as population index.");

public:
    using Base = CompartmentalModel<Comp, Pop, Params>;
    /**
     * @brief default constructor, forwarding args to Base constructor
     */
    template <class... Args>
    CompartmentalModel<Comp, Pop, Params, FlowChart<Flows...>>(Args... args)
        : CompartmentalModel<Comp, Pop, Params>(args...)
        , m_comp_factor(comp_factor())
        , m_flow_values((this->populations.numel() / static_cast<size_t>(Comp::Count)) * FlowChart<Flows...>().size())
    {
    }

    // Note: use get_flow_index when accessing flows
    virtual void get_flows(Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<const Eigen::VectorXd> /*y*/, double /*t*/,
                           Eigen::Ref<Eigen::VectorXd> /*flows*/) const {};

    /**
     * @brief compute the right-hand-side of the ODE dydt = f(y, t) from flow values
     *
     * This function is generated at compile time depending on the FlowChart tparam of the model.
     *
     * @param flows the current flow values (as calculated by get_flows) as a flat array
     * @param dydt a reference to the calculated output
     */
    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> flows, Eigen::Ref<Eigen::VectorXd> dydt) const
    {
        dydt.setZero();
        // create a new population index, by copying its size
        // Note: we would like to have Index{0, ..., 0}, but since there is no default ctor, this is not easy to get.
        //       -> all entries of index must be overwritten!
        auto index = this->populations.size();
        get_rhs_impl(flows, dydt, index);
    }

    /**
     * @brief compute the right-hand-side of the ODE dydt = f(y, t)
     *
     * This function is uses get_flow(..., flows) and get_derivatives(flows, dydt) to provide the
     * same interface as a CompartmentalModel without flows
     *
     * @param pop the current population of the model as a flat array
     * @param y the current state of the model as a flat array
     * @param t the current time
     * @param dydt a reference to the calculated output
     */
    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override final
    {
        m_flow_values.setZero();
        get_flows(pop, y, t, m_flow_values);
        get_derivatives(m_flow_values, dydt);
    }

    /**
     * @brief initial values for flows.
     * This can be used as initial conditions in an ODE solver. By default, this is a zero vector
     * @return the initial flows
     */
    Eigen::VectorXd get_initial_flows() const
    {
        return Eigen::VectorXd(FlowChart<Flows...>().size());
    }

    /**
     * @brief a flat index into an array of flows (as computed by get_flows), given the indices of each category
     * @param indices the custom indices for each category
     * @tparam Source the source of a flow
     * @tparam Target the target of a flow
     * @return a flat index into a data structure storing flow values
     */
    template <Comp Source, Comp Target>
    size_t get_flow_index(const FlowIndex& indices) const
    {
        // Calculation of flat indices (as in this->populations.get_flat_index) use the formula
        // flat_index = ((((I_{1}) * D_{2} + I_{2}) * D_{3} + I_{3}) ... ) * D_{n} + I_{n}
        //            = I_{n} + I_{n-1} * D_{n} + I_{n-2} * D_{n} * D_{n-1} + ... + I_{1} * D_{n} * ... * D_{2}
        // for indices (I_{1}, ... , I_{n}) and dimension (D_{1}, ... , D_{n}).
        //
        // The flows use their position in the FlowChart instead of Comp. As this dimension (the number of flows) may
        // be larger or smaller (e.g. a S->I->R model has 3 Comp's, but 2 Flows) than Comp::Count, we cannot simply use
        // this->populations.get_flat_index.
        //
        // Instead, we omit Comp from PopIndex (getting FlowIndex), calculate an without the Comp index or Dimension,
        // i.e. we omit I_j and D_j corresponding to Comp from the formula above (via flatten_index_by_tags).
        // Then, we add the flow position I_{flow} via flow_index = flat_index * D_{flow} + I_{flow}, where
        // D_{flow} is the number of flows.

        return flatten_index_by_tags<FlowIndex>(indices, this->populations.size()) * FlowChart<Flows...>().size() +
               FlowChart<Flows...>().template get<Flow<Comp, Source, Target>>();
    }

    /**
     * @brief a flat index into an array of flows (as computed by get_flows), if the only used category in Pop is Comp
     * @tparam Source the source of a flow
     * @tparam Target the target of a flow
     * @return a flat index into a data structure storing flow values
     */
    template <Comp Source, Comp Target>
    constexpr size_t get_flow_index() const
    {
        static_assert(std::is_same<FlowIndex, Index<>>::value, "Other indizes must be specified");
        return FlowChart<Flows...>().template get<Flow<Comp, Source, Target>>();
    }

private:
    size_t m_comp_factor;
    mutable Eigen::VectorXd m_flow_values;

    // compute only the factor used with "Comp" in this->population.get_flat_index
    template <size_t I = PopIndex::size - 1>
    inline std::enable_if_t<(I > details::IndexPosition<Comp, PopIndex>::value), size_t> comp_factor()
    {
        static_assert(I > 0, "Population must contain Compartments Category.");
        // see get_flow_index for more information on flat indizes
        return static_cast<size_t>(get<I>(this->populations.size())) * comp_factor<I - 1>();
    }

    // compute only the factor used with "Comp" in this->population.get_flat_index
    template <size_t I = PopIndex::size - 1>
    inline std::enable_if_t<(I <= details::IndexPosition<Comp, PopIndex>::value), size_t> comp_factor()
    {
        return 1;
    }

    // recursively compute the derivatives/rhs of the model equations.
    // Uses runtime loops for non-"Comp" indices via the get_rhs_outer_impl, and compile-time recursion for Flows
    // (which provide "Comp" indices as their source/target) via the get_rhs_inner_impl, called from "outer"
    template <template <class...> class Index, class... Categories>
    void get_rhs_impl(Eigen::Ref<const Eigen::VectorXd> flows, Eigen::Ref<Eigen::VectorXd> rhs,
                      Index<Categories...>& I) const
    {
        get_rhs_outer_impl<void, Categories...>(flows, rhs, I);
    }

    // compute rhs using flows - recursive function that loops over each Category in Populations, skipping Comp
    // this function handles general Categories, e.g. AgeGroups
    template <class Dummy, class Category, class... Categories>
    inline std::enable_if_t<!std::is_same<Category, Comp>::value>
    get_rhs_outer_impl(Eigen::Ref<const Eigen::VectorXd> flows, Eigen::Ref<Eigen::VectorXd> rhs, PopIndex& I) const
    {
        const size_t size = static_cast<size_t>(get<Category>(this->populations.size()));
        // loop over the given Category, setting the corresponding entry in I accordingly
        for (size_t i = 0; i < size; ++i) {
            get<Category>(I) = Index<Category>(i);
            get_rhs_outer_impl<Dummy, Categories...>(flows, rhs, I);
        }
    }

    // compute rhs using flows - recursive function that loops over each Category in Populations, skipping Comp
    // this function handles Category==Comp (by not doing anything. get_rhs_inner_impl will deal with it)
    template <class Dummy, class Category, class... Categories>
    inline std::enable_if_t<std::is_same<Category, Comp>::value>
    get_rhs_outer_impl(Eigen::Ref<const Eigen::VectorXd> flows, Eigen::Ref<Eigen::VectorXd> rhs, PopIndex& I) const
    {
        get_rhs_outer_impl<Dummy, Categories...>(flows, rhs, I);
    }

    // compute rhs using flows - recursive function that loops over each Category in Populations, skipping Comp
    // this function pre-calculates flat indices, calls get_rhs_inner_impl, and stops the recursion
    template <class Dummy>
    inline void get_rhs_outer_impl(Eigen::Ref<const Eigen::VectorXd> flows, Eigen::Ref<Eigen::VectorXd> rhs,
                                   PopIndex& I) const
    {
        // set Comp = 0, so that the inner recursion can add their relevant index, using m_comp_factor
        get<Comp>(I) = Index<Comp>(0);
        // perform inner recursion, i.e. calculate the rhs for all flows, given index I
        get_rhs_inner_impl<0>(flows, rhs, this->populations.get_flat_index(I),
                              flatten_index_by_tags<FlowIndex>(I, this->populations.size()));
        // end recursion
    }

    // recursively compute the derivatives/rhs of the model equations for a given Flow and non-"Comp" (flat) indices
    template <size_t flow_id>
    inline std::enable_if_t<flow_id != sizeof...(Flows)>
    get_rhs_inner_impl(Eigen::Ref<const Eigen::VectorXd> flows, Eigen::Ref<Eigen::VectorXd> rhs, const size_t& flat_pop,
                       const size_t& flat_flow) const
    {
        // see get_flow_index for more info on these index calculations
        constexpr size_t source = static_cast<size_t>(FlowChart<Flows...>().template get<flow_id>().source);
        constexpr size_t target = static_cast<size_t>(FlowChart<Flows...>().template get<flow_id>().target);
        // flat_flow and flat_pop already contain index information for all non-Comp Categories, so
        // we only need to add flow_id and source/target respectively (with the correct factors)
        const size_t flow_index = flow_id + flat_flow * FlowChart<Flows...>().size();
        // flat_flow already contains indices for all non-Comp Categories, so we only need to add the flow_id here
        rhs[flat_pop + source * m_comp_factor] -= flows[flow_index];
        rhs[flat_pop + target * m_comp_factor] += flows[flow_index];
        get_rhs_inner_impl<flow_id + 1>(flows, rhs, flat_pop, flat_flow);
    }

    template <size_t flow_id>
    inline std::enable_if_t<flow_id == sizeof...(Flows)> get_rhs_inner_impl(Eigen::Ref<const Eigen::VectorXd>,
                                                                            Eigen::Ref<Eigen::VectorXd>, const size_t&,
                                                                            const size_t&) const
    {
        // end recursion
    }
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
