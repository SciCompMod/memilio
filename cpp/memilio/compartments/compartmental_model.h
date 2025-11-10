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
#ifndef MIO_COMPARTMENTS_COMPARTMENTAL_MODEL_H
#define MIO_COMPARTMENTS_COMPARTMENTAL_MODEL_H

#include "memilio/config.h"
#include "memilio/math/eigen.h"
#include <concepts>

namespace mio
{

template <class T>
concept HasCheckConstraints = requires(T t) {
    { t.check_constraints() } -> std::same_as<bool>;
};

template <class T>
concept HasApplyConstraints = requires(T t) {
    { t.apply_constraints() } -> std::same_as<bool>;
};

/**
 * @brief CompartmentalModel is a template for a compartmental model for an
 * array of initial populations and a parameter set.
 * @tparam FP A floating point type, e.g., double.
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

    CompartmentalModel(const CompartmentalModel&)            = default;
    CompartmentalModel(CompartmentalModel&&)                 = default;
    CompartmentalModel& operator=(const CompartmentalModel&) = default;
    CompartmentalModel& operator=(CompartmentalModel&&)      = default;
    virtual ~CompartmentalModel()                            = default;

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
        if constexpr (HasApplyConstraints<ParameterSet>) {
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
        if constexpr (HasCheckConstraints<ParameterSet>) {
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
 * @brief Concept to check if a type is a valid compartment model.
 * Note that Model must be the first template argument 
 * @tparam Model a type that may or may not be a compartment model.
 * @tparam FP, floating point type, e.g., double.
 */
template <class Model, typename FP>
concept IsCompartmentalModel =
    requires(Model m, Eigen::Ref<const Eigen::VectorX<FP>> pop_y, Eigen::Ref<Eigen::VectorX<FP>> dydt, FP t) {
        { m.get_initial_values() } -> std::convertible_to<Eigen::VectorX<FP>>;
        m.eval_right_hand_side(pop_y, pop_y, t, dydt);
    };

} // namespace mio

#endif // MIO_COMPARTMENTS_COMPARTMENTAL_MODEL_H
