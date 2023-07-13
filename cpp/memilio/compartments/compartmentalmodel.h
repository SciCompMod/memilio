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
#ifndef MIO_COMPARTMENTALMODEL_H
#define MIO_COMPARTMENTALMODEL_H

#include "memilio/config.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/type_chart.h"
#include "memilio/utils/flow.h"
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

namespace details
{
// The following functions are not defined anywhere. Their use is to provide type conversions via their return type.

// Function declaration used to remove OmittedTag from the type list of a tuple.
// First a list of tuples is generated for each Tag in Tags, where the tuple is either of type tuple<Tag>, or if
// Tag == OmittedTag, of type tuple<>. This list is then concatonated, effectively removing OmittedTag.
template <class OmittedTag, class... Tags>
decltype(std::tuple_cat(std::declval<typename std::conditional<std::is_same<OmittedTag, Tags>::value, std::tuple<>,
                                                               std::tuple<Tags>>::type>()...))
    filter_tuple(std::tuple<Tags...>);

// Function declaration used to replace type T by std::tuple.
template <template <class...> class T, class... Args>
std::tuple<Args...> as_tuple(T<Args...>);

// Function declaration used to replace std::tuple by type T.
template <template <class...> class T, class... Args>
T<Args...> as_index(std::tuple<Args...>);

// Remove all occurrences of OmittedTag from the types in a std::tuple<types...>.
template <class OmittedTag, class Tuple>
using filtered_tuple_t = decltype(filter_tuple<OmittedTag>(std::declval<Tuple>()));

// Remove all occurrences of OmittedTag from the types in an Index = IndexTemplate<types...>.
template <class OmittedTag, template <class...> class IndexTemplate, class Index>
using filtered_index_t = decltype(as_index<IndexTemplate>(
    std::declval<filtered_tuple_t<OmittedTag, decltype(as_tuple(std::declval<Index>()))>>()));

} //namespace details

/**
 * @brief CompartmentalModel with flows. 
 *
 * Inherits from CompartmentalModel extending it with a TypeChart (the optional template parameter). 
 * @see CompartmentalModel
 *
 * Instead of directly computing the right-hand-side with get_derivatives, first the flows defined by TypeChart are
 * computed, and then added to/subtracted from their source/target compartment.
 */
template <class Comp, class Pop, class Params, class... Flows>
struct CompartmentalModel<Comp, Pop, Params, TypeChart<Flows...>> : public CompartmentalModel<Comp, Pop, Params> {
    using PopIndex = typename Pop::Index;
    // FlowIndex is the same as PopIndex without the category Comp. It is used as argument type for get_flow_index,
    // since a flow is used to determine only the compartment (i.e. Index<Comp>) for the population of the model.
    // The remaining indices (those contained in FlowIndex) must still be provided to get an entry of population.
    // This approach only works with exactly one category of type Comp in PopIndex (hence the assertion below).
    using FlowIndex = details::filtered_index_t<Comp, Index, PopIndex>;
    // Enforce that Comp is a unique Category of PopIndex, since we use Flows (via their source/target) to provide
    // Comp indices for the population.
    static_assert(FlowIndex::size == PopIndex::size - 1, "Compartments must be used exactly once as population index.");

public:
    using Base      = CompartmentalModel<Comp, Pop, Params>;
    using FlowChart = TypeChart<Flows...>;
    /**
     * @brief Default constructor, forwarding args to Base constructor.
     */
    template <class... Args>
    CompartmentalModel<Comp, Pop, Params, FlowChart>(Args... args)
        : CompartmentalModel<Comp, Pop, Params>(args...)
        , m_comp_factor(comp_factor())
        , m_flow_index_dimensions(reduce_index(this->populations.size()))
        , m_flow_values((this->populations.numel() / static_cast<size_t>(Comp::Count)) * FlowChart().size())
    {
    }

    // Note: use get_flow_index when accessing flows
    // Note: by convention, we compute incoming flows, thus entries in flows must be non-negative
    virtual void get_flows(Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<const Eigen::VectorXd> /*y*/, double /*t*/,
                           Eigen::Ref<Eigen::VectorXd> /*flows*/) const {};

    /**
     * @brief Compute the right-hand-side of the ODE dydt = f(y, t) from flow values.
     *
     * This function is generated at compile time depending on the template parameters Flows and Pop of the model.
     *
     * @param[in] flows The current flow values (as calculated by get_flows) as a flat array.
     * @param[out] dydt A reference to the calculated output.
     */
    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> flows, Eigen::Ref<Eigen::VectorXd> dydt) const
    {
        dydt.setZero();
        for (size_t i = 0; i < this->populations.numel() / static_cast<size_t>(Comp::Count); i++) {
            // see get_flow_index for more info on these index calculations.
            // the Comp-index is added in get_rhs_impl; the flat PopIndex is
            //   flat_pop_index = (i % m_comp_factor) + (i / m_comp_factor) * m_comp_factor * Comp::Count;
            // mind the integer division.
            get_rhs_impl<0>(flows, dydt,
                            i + (i / m_comp_factor) * m_comp_factor * (static_cast<size_t>(Comp::Count) - 1),
                            i * FlowChart().size());
        }
    }

    /**
     * @brief Compute the right-hand-side f(y, t) of the ODE and store it in dydt.
     *
     * This function uses get_flow(..., flows) and get_derivatives(flows, dydt) to provide the
     * same interface as a CompartmentalModel without flows.
     *
     * @param[in] pop The current population of the model as a flat array.
     * @param[in] y The current state of the model as a flat array.
     * @param[in] t The current time.
     * @param[out] dydt A reference to the calculated output.
     */
    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override final
    {
        m_flow_values.setZero();
        get_flows(pop, y, t, m_flow_values);
        get_derivatives(m_flow_values, dydt);
    }

    /**
     * @brief Initial values for flows.
     * This can be used as initial conditions in an ODE solver. By default, this is a zero vector.
     * @return The initial flows.
     */
    Eigen::VectorXd get_initial_flows() const
    {
        return Eigen::VectorXd::Zero((this->populations.numel() / static_cast<size_t>(Comp::Count)) *
                                     FlowChart().size());
    }

    /**
     * @brief S flat index into an array of flows (as computed by get_flows), given the indices of each category.
     *
     * Calculation of flat indices (as in this->populations.get_flat_index) use the formula
     * flat_index = ((((I_{1}) * D_{2} + I_{2}) * D_{3} + I_{3}) ... ) * D_{n} + I_{n}
     * = I_{n} + I_{n-1} * D_{n} + I_{n-2} * D_{n} * D_{n-1} + ... + I_{1} * D_{n} * ... * D_{2}
     * for indices (I_{1}, ... , I_{n}) and dimension (D_{1}, ... , D_{n}).
     *   
     * The flows use their position in the FlowChart instead of the template argument Comp from the base class,
     * i.e. they use the Population::Index type without Comp. The position is determined by the Source and Target 
     * template arguments instead. As this dimension (the number of flows) may be larger or smaller
     * (e.g. a S->I->R model has 3 Comp's, but 2 Flows) than Comp::Count, we cannot simply use
     * this->populations.get_flat_index.
     *   
     * Instead, we omit Comp from PopIndex (getting FlowIndex), calculate the flat_index without the Comp index or
     * Dimension, i.e. we omit I_j and D_j corresponding to Comp from the formula above. We do this by calling
     * flatten_index with a FlowIndex, with dimensions derived from Pop via reduce_index.
     * Then, we add the flow position I_{flow} via flow_index = flat_index * D_{flow} + I_{flow}, where
     * D_{flow} is the number of flows.
     *
     * @param[in] indices The custom indices for each category.
     * @tparam Source the source of a flow.
     * @tparam Target the target of a flow.
     * @return A flat index into a data structure storing flow values.
     */
    template <Comp Source, Comp Target>
    size_t get_flow_index(const FlowIndex& indices) const
    {
        return flatten_index(indices, m_flow_index_dimensions) * FlowChart().size() +
               FlowChart().template get<Flow<Comp, Source, Target>>();
    }

    /**
     * @brief A flat index into an array of flows (as computed by get_flows), if the only used category in Pop is Comp.
     * @tparam Source the source of a flow.
     * @tparam Target the target of a flow.
     * @return A flat index into a data structure storing flow values.
     */
    template <Comp Source, Comp Target>
    constexpr size_t get_flow_index() const
    {
        static_assert(std::is_same<FlowIndex, Index<>>::value, "Other indices must be specified");
        return FlowChart().template get<Flow<Comp, Source, Target>>();
    }

private:
    size_t m_comp_factor; ///< The factor used for the category Comp in this->population.get_flat_index.
    FlowIndex m_flow_index_dimensions; ///< The dimensions of a FlowIndex.
    mutable Eigen::VectorXd m_flow_values; ///< Cache to avoid allocation in get_derivatives (using get_flows).

    /**
     * @brief Reduce a PopIndex to a FlowIndex.
     * This function should only be called by `FlowIndex reduce_index(const PopIndex&)`.
     */
    template <class... Categories>
    FlowIndex reduce_index(const PopIndex& i, mio::Tag<mio::Index<Categories...>>) const
    {
        return Index<Categories...>{get<Categories>(i)...};
    }

    /**
     * @brief Reduce a PopIndex to a FlowIndex.
     * @param[in] i Any PopIndex.
     * @return The FlowIndex gained from i by omitting the compartment.
     */
    FlowIndex reduce_index(const PopIndex& i) const
    {
        // since FlowIndex may not be trivially constructible, we use mio::Tag<FlowIndex> to determine its Categories
        return reduce_index(i, mio::Tag<FlowIndex>());
    }

    /**
     * @brief Recursively compute the factor used with the Category "Comp" in this->population.get_flat_index.
     * @tparam I Recursion index, do not change default value.
     * @return Factor used to add a compartment index to a flat index.
     */
    template <size_t I = PopIndex::size - 1>
    inline std::enable_if_t<(I > details::IndexPosition<Comp, PopIndex>::value), size_t> comp_factor()
    {
        static_assert(I > 0, "Population must contain Compartments Category.");
        // see get_flow_index for more information on flat indizes
        return static_cast<size_t>(get<I>(this->populations.size())) * comp_factor<I - 1>();
    }

    /// @brief Ends recursion of `size_t comp_factor()`.
    template <size_t I = PopIndex::size - 1>
    inline std::enable_if_t<(I <= details::IndexPosition<Comp, PopIndex>::value), size_t> comp_factor()
    {
        return 1;
    }

    /**
     * @brief Compute the derivatives of the compartments.
     * Compute the derivatives/rhs of the model equations for a given flow and non-"Comp" (flat) indices.
     * Uses recursion, which is resolved at compile time depending on the FlowChart.
     * @param[in] flows Current change in flows, as computed by get_flows.
     * @param[out] rhs The derivatives of the model equations.
     * @param[in] flat_pop Flattened PopIndex without contribution from the compartment index.
     * @param[in] flat_flow Flattened FlowIndex.
     * @tparam FlowID The index of a flow in FlowChart.
     */
    template <size_t FlowID>
    inline std::enable_if_t<FlowID != sizeof...(Flows)>
    get_rhs_impl(Eigen::Ref<const Eigen::VectorXd> flows, Eigen::Ref<Eigen::VectorXd> rhs, const size_t& flat_pop,
                 const size_t& flat_flow) const
    {
        constexpr size_t source = static_cast<size_t>(FlowChart().template get<FlowID>().source);
        constexpr size_t target = static_cast<size_t>(FlowChart().template get<FlowID>().target);
        // see get_flow_index for more info on these index calculations.
        // flat_flow and flat_pop already contain index information for all non-Comp Categories, so
        // we only need to add FlowID and source/target respectively (with the correct factors).
        rhs[flat_pop + source * m_comp_factor] -= flows[flat_flow + FlowID];
        rhs[flat_pop + target * m_comp_factor] += flows[flat_flow + FlowID];
        get_rhs_impl<FlowID + 1>(flows, rhs, flat_pop, flat_flow);
    }

    /// @brief Ends recursion of `void get_rhs_impl(...)`.
    template <size_t FlowID>
    inline std::enable_if_t<FlowID == sizeof...(Flows)>
    get_rhs_impl(Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<Eigen::VectorXd>, const size_t&, const size_t&) const
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

#endif // MIO_COMPARTMENTALMODEL_H
