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

    CompartmentalModel(const CompartmentalModel&) = default;
    CompartmentalModel(CompartmentalModel&&)      = default;
    CompartmentalModel& operator=(const CompartmentalModel&) = default;
    CompartmentalModel& operator=(CompartmentalModel&&) = default;
    virtual ~CompartmentalModel()                       = default;

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
 * Inherited class from CompartmentalModel extended with flows. 
 * The flows are first computed and then assigned to their source and target,
 * allowing us to determine the population in the infection states.
 */
template <class Comp, class Pop, class Params, class... Flows>
struct CompartmentalModel<Comp, Pop, Params, FlowChart<Flows...>> : public CompartmentalModel<Comp, Pop, Params> {
    using FlowIndex = details::filtered_index_t<Comp, Index, typename Pop::Index>;
    static_assert(FlowIndex::size == Pop::Index::size - 1,
                  "Compartments must be used exactly once as population index.");

public:
    template <class... Args>
    CompartmentalModel<Comp, Pop, Params, FlowChart<Flows...>>(Args... args)
        : CompartmentalModel<Comp, Pop, Params>(args...)
        , m_flow_dims(make_dims(this->populations))
        , m_flow_flat_dim(product(m_flow_dims))
        , m_flow_values(m_flow_flat_dim * FlowChart<Flows...>().size())
    {
    }

    virtual void get_flows(Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<const Eigen::VectorXd> /*y*/, double /*t*/,
                           Eigen::Ref<Eigen::VectorXd> /*dydt*/) const {};

    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> flows, Eigen::Ref<Eigen::VectorXd> dydt) const
    {
        dydt.setZero();
        const size_t num_compartments = dydt.size() / m_flow_flat_dim;
        for (size_t i = 0; i < m_flow_flat_dim; i++) {
            get_rhs_impl<0>(flows, dydt, i, num_compartments);
        }
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override final
    {
        m_flow_values.setZero();
        get_flows(pop, y, t, m_flow_values);
        get_derivatives(m_flow_values, dydt);
    }

    // returns initial conditions when simulation with flows. Here, we return an empty eigen vector with the same size as number of flows.
    Eigen::VectorXd get_initial_flows() const
    {
        return Eigen::VectorXd(FlowChart<Flows...>().size());
    }

protected:
    template <Comp Source, Comp Target>
    size_t get_flow_index(const FlowIndex& f) const
    {
        return flatten_index(f, m_flow_dims) +
               m_flow_flat_dim * FlowChart<Flows...>().template get<Flow<Comp, Source, Target>>();
    }

    template <Comp Source, Comp Target>
    constexpr size_t get_flow_index() const
    {
        static_assert(std::is_same<FlowIndex, Index<>>::value, "Other indizes must be specified");
        return FlowChart<Flows...>().template get<Flow<Comp, Source, Target>>();
    }

private:
    const FlowIndex m_flow_dims;
    const size_t m_flow_flat_dim;
    mutable Eigen::VectorXd m_flow_values;

    template <class Index, class Tag, class... Tags>
    constexpr std::enable_if_t<!std::is_same<Tag, Comp>::value, FlowIndex>
    make_flow_index_impl(const Index& other) const
    {
        auto i      = make_flow_index_impl<Index, Tags...>(other);
        get<Tag>(i) = other.template size<Tag>();
        return i;
    }

    template <class Index, class Tag, class... Tags>
    constexpr std::enable_if_t<std::is_same<Tag, Comp>::value, FlowIndex> make_flow_index_impl(const Index& other) const
    {
        return make_flow_index_impl<Index, Tags...>(other);
    }

    template <class Index>
    constexpr FlowIndex make_flow_index_impl(const Index&) const
    {
        return FlowIndex();
    }

    template <class Type, class... Tags>
    constexpr FlowIndex make_dims(CustomIndexArray<Type, Tags...> i) const
    {
        return make_flow_index_impl<CustomIndexArray<Type, Tags...>, Tags...>(i);
    }

    // recursive implemention of mapping  flow source/target to the corresponding rhs entry.
    template <size_t i>
    inline std::enable_if_t<i != sizeof...(Flows)>
    get_rhs_impl(Eigen::Ref<const Eigen::VectorXd> flows, Eigen::Ref<Eigen::VectorXd> rhs, const size_t index_age_group,
                 const size_t num_compartments) const
    {
        size_t indx_source =
            static_cast<size_t>(FlowChart<Flows...>().template get<i>().source) + num_compartments * index_age_group;
        size_t indx_target =
            static_cast<size_t>(FlowChart<Flows...>().template get<i>().target) + num_compartments * index_age_group;

        rhs[indx_source] -= flows[index_age_group + m_flow_flat_dim * i];
        rhs[indx_target] += flows[index_age_group + m_flow_flat_dim * i];
        get_rhs_impl<i + 1>(flows, rhs, index_age_group, num_compartments);
    }

    template <size_t i>
    inline std::enable_if_t<i == sizeof...(Flows)>
    get_rhs_impl(Eigen::Ref<const Eigen::VectorXd>, Eigen::Ref<Eigen::VectorXd>, const size_t, const size_t) const
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
