/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Anna Wendler, Lena Ploetzke
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
#ifndef STATEAGEFUNCTION_H
#define STATEAGEFUNCTION_H

#include "memilio/config.h"
#include "memilio/utils/parameter_set.h"
#include "ide_secir/infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/math/smoother.h"
#include "memilio/math/floating_point.h"
#include "memilio/epidemiology/uncertain_matrix.h"

namespace mio
{
/**************************
* Define StateAgeFunction *
***************************/

/**
 * @brief A generic function depending on the state age, i.e. the time already spent in some #InfectionState.
 * This is an abstract class and cannot be instantiated on its own.
 *
 * By construction, StateAgeFunction%s S(x) are only defined for non-negative values x>=0 with S(0)>=0. In order to
 * approximate the support (size) of a StateAgeFunction it is thus sufficient to evaluate the function for positive x.
 *
 * Derived StateAgeFunctions can be used for two types of functionality. 
 *  a) Monotonously decreasing functions describing the share of individuals that have not yet transitioned to the next
 *     #InfectionState. These functions are denoted as #TransitionDistributions since 1 - TransitionDistribution 
 *     represents a cumulative distribution function.
 *  b) Arbitrary non-negative functions used for parameters such as TransmissionProbabilityOnContact.
 * 
 * Derived classes must implement the eval method which implements the actual function that is evaluated at some state age.
 * This function can depend on one parameter.
 *
 * The derived classes must also implement the clone_impl method which allows to deepcopy the derived class.
 * 
 * The get_support_max method is virtual and implements a basic version to determine the maximum of the support. 
 * For some derived classes there is a more efficient way (see e.g., SmootherCosine) to do this which is 
 * why it can be overridden. The base class implementation uses the fact that the StateAgeFunction is monotonously 
 * decreasing. This is no limitation as the support is only needed for StateAgeFunctions of Type a) as given above.
 * For classes of type b) a dummy implementation logging an error and returning -2 for get_support_max() should be implemented.
 *
 * See ExponentialDecay, SmootherCosine and ConstantFunction for examples of derived classes.
 */
struct StateAgeFunction {

    /**
     * @brief Constructs a new StateAgeFunction object
     * 
     * @param[in] init_parameter Specifies the initial function parameter of the function.
     */
    StateAgeFunction(ScalarType init_parameter)
        : m_parameter{init_parameter}
        , m_support_max{-1.} // initialize support maximum as not set
        , m_support_tol{-1.} // initialize support tolerance as not set
    {
    }

    /**
     * @brief Virtual destructor.
     */
    virtual ~StateAgeFunction() = default;

    /**
     * @brief Copy constructor.
     */
    StateAgeFunction(const StateAgeFunction& other) = default;

    /**
     * @brief Move constructor.
     */
    StateAgeFunction(StateAgeFunction&& other) = default;

    /**
     * @brief Copy assignment operator.
     */
    StateAgeFunction& operator=(const StateAgeFunction& other) = default;

    /**
     * @brief Move assignment operator.
     */
    StateAgeFunction& operator=(StateAgeFunction&& other) = default;

    /**
     * @brief Comparison operator.
     */
    bool operator==(const StateAgeFunction& other) const
    {
        return (typeid(*this).name() == typeid(other).name() && m_parameter == other.get_parameter());
    }

    /**
     * @brief Here a pure virtual function is defined that depends on the state_age. 
     *
     * The defined function ususally depends on some function parameter.
     * 
     * @param[in] state_age Time at which the function is evaluated.
     */
    virtual ScalarType eval(ScalarType state_age) = 0;

    /**
     * @brief Get the m_parameter object
     * 
     * Can be used to access the m_parameter object, which specifies the used function.
     * 
     * @return ScalarType 
     */
    ScalarType get_parameter() const
    {
        return m_parameter;
    }

    /**
     * @brief Set the m_parameter object.
     * 
     * Can be used to set the m_parameter object, which specifies the used function.
     * The maximum support of a function may be costly to evaluate. In order to not always reevaluate or recompute the
     * support when the user asks for it, a cached value is used. If m_support_max is set to -1, the cached value is
     * deleted and a recomputation is done the next time the user asks for the support. As the support (potentially)
     * depends on the m_parameter object, the cached value has to be deleted. For details see get_support_max().
     *
     *@param[in] new_parameter New parameter for StateAgeFunction.
     */
    void set_parameter(ScalarType new_parameter)
    {
        m_parameter = new_parameter;

        m_support_max = -1.;
    }

    /**
     * @brief Computes the maximum of the support of the function using the time step size dt and some tolerance tol. 
     * 
     * This is a basic version to determine the maximum of the support of a monotonously decreasing function,
     * evaluating the function at each time point until it reaches the boundaries of the support.
     *
     * For some specific derivations of StateAgeFunction%s there are more efficient ways to determine the 
     * support which is why this member function is virtual and can be overridden (see, e.g., SmootherCosine).
     * The maximum of the support is only needed for StateAgeFunction%s that are used as TransitionDistribution%s. 
     *
     * @param[in] dt Time step size at which function will be evaluated. 
     * @param[in] tol Tolerance used for cutting the support if the function value falls below. 
     * @return ScalarType support_max
     */
    virtual ScalarType get_support_max(ScalarType dt, ScalarType tol = 1e-10)
    {
        ScalarType support_max = 0.;

        if (!floating_point_equal(m_support_tol, tol, 1e-14) || floating_point_equal(m_support_max, -1., 1e-14)) {
            while (eval(support_max) >= tol) {
                support_max += dt;
            }

            m_support_max = support_max;
            m_support_tol = tol;
        }

        return m_support_max;
    }

    /**
     * @brief Get type of StateAgeFunction, i.e.which derived class is used.
     * 
     * @param[out] string 
     */
    std::string get_state_age_function_type() const
    {
        return typeid(*this).name();
    }

    /**
     * @brief Clones unique pointer to a StateAgeFunction.
     * 
     * Calls the clone_impl method that is implemented by every derived class.
     * 
     * @return std::unique_ptr<StateAgeFunction> unique pointer to a StateAgeFunction
     */
    std::unique_ptr<StateAgeFunction> clone() const
    {
        return std::unique_ptr<StateAgeFunction>(clone_impl());
    }

protected:
    /**
     * @brief Pure virtual method that implements cloning.
     */
    virtual StateAgeFunction* clone_impl() const = 0;

    ScalarType m_parameter; ///< Parameter for function in derived class.
    ScalarType m_support_max; ///< Maximum of the support of the function.
    ScalarType m_support_tol; ///< Tolerance for computation of the support.
};

/**************************************
* Derived classes of StateAgeFunction *
***************************************/

/**
 * @brief Class that defines an exponential decay function depending on the state age.
 */
struct ExponentialDecay : public StateAgeFunction {

    /**
     * @brief Constructs a new ExponentialDecay object
     * 
     * @param[in] init_parameter Specifies the initial function parameter of the function.
     */
    ExponentialDecay(ScalarType init_parameter)
        : StateAgeFunction(init_parameter)
    {
    }

    /**
     * @brief Defines exponential decay function depending on state_age.
     *
     * m_parameter defines how fast the exponential function decays.
     * 
     * @param[in] state_age Time at which the function is evaluated.
     * @return Evaluation of the function at state_age. 
     */
    ScalarType eval(ScalarType state_age) override
    {
        return std::exp(-m_parameter * state_age);
    }

protected:
    /**
     * @brief Implements clone for ExponentialDecay.
     * 
     * @return Pointer to StateAgeFunction.
     */
    StateAgeFunction* clone_impl() const override
    {
        return new ExponentialDecay(*this);
    }
};

/**
 * @brief Class that defines an smoother_cosine function depending on the state age.
 */
struct SmootherCosine : public StateAgeFunction {

    /**
     * @brief Constructs a new SmootherCosine object
     * 
     * @param[in] init_parameter specifies the initial parameter of the function.
     */
    SmootherCosine(ScalarType init_parameter)
        : StateAgeFunction(init_parameter)
    {
    }

    /**
     * @brief Defines smoother cosine function depending on state_age.
     *
     * Used function goes through points (0,1) and (m_parameter,0) and is interpolated in between using a smoothed cosine function.
     * 
     * @param[in] state_age Time at which the function is evaluated.
     * @return Evaluation of the function at state_age. 
     */
    ScalarType eval(ScalarType state_age) override
    {
        return smoother_cosine(state_age, 0.0, m_parameter, 1.0, 0.0);
    }

    /**
     * @brief Computes the maximum of the support of the function. 
     * 
     * For SmootherCosine, the maximum of the support is equal to the function parameter.
     *
     * @param[in] dt Time step size. 
     * @param[in] tol Tolerance used for cutting the support if the function value falls below. 
     * @return ScalarType support_max
     */
    ScalarType get_support_max(ScalarType dt, ScalarType tol = 1e-10) override
    {
        unused(dt);
        unused(tol);
        m_support_max = m_parameter;
        return m_support_max;
    }

protected:
    /**
     * @brief Clones unique pointer to a StateAgeFunction.
     * 
     * @return std::unique_ptr<StateAgeFunction> unique pointer to a StateAgeFunction
     */
    StateAgeFunction* clone_impl() const override
    {
        return new SmootherCosine(*this);
    }
};

/**
 * @brief Class that defines a constant function.
 */
struct ConstantFunction : public StateAgeFunction {

    /**
     * @brief Constructs a new ConstantFunction object
     * 
     * @param init_parameter specifies value of the constant function.
     */
    ConstantFunction(ScalarType init_parameter)
        : StateAgeFunction(init_parameter)
    {
    }

    /**
     * @brief Defines constant function.
     *
     *  The function parameter defines the value of the function. 
     *
     * @param state_age Time at which the function is evaluated.
     * @return Evaluation of the function at state_age. 
     */
    ScalarType eval(ScalarType state_age) override
    {
        unused(state_age);
        return m_parameter;
    }

    /**
     * @brief Computes the maximum of the support of the function. 
     * 
     * For ConstantFunction the maximum of the support would be infinity. This is why we do not want to use it 
     * as a TransitionDistribution and getting the maximum of the support doe not make sense. 
     *
     * @param[in] dt Time step size. 
     * @param[in] tol Tolerance used for cutting the support if the function value falls below. 
     * @return ScalarType support_max
     */
    ScalarType get_support_max(ScalarType dt, ScalarType tol = 1e-10) override
    {
        // In case of a ConstantFunction we would have support_max = infinity
        // This type of function is not suited to be a TransitionDistribution
        // Log error and return -2.

        unused(dt);
        unused(tol);
        m_support_max = -2.;

        log_error("This function is not suited to be a TransitionDistribution. Do not call in case of StateAgeFunctions"
                  "of type b); see documentation of StateAgeFunction Base class.");

        return m_support_max;
    }

protected:
    /**
     * @brief Clones unique pointer to a StateAgeFunction.
     * 
     * @return std::unique_ptr<StateAgeFunction> unique pointer to a StateAgeFunction
     */
    StateAgeFunction* clone_impl() const override
    {
        return new ConstantFunction(*this);
    }
};

/*********************************
* Define StateAgeFunctionWrapper *
*********************************/

/**
 * @brief Wrapper around StateAgeFunction so that one can work with an arbitrary StateAgeFunction.
 *
 * This way we can define e.g. the parameter TransmissionProbabilityOnContact as type StateAgeFunctionWrapper
 * and set it with a specific StateAgeFunction for each example.  
 * 
 * Example from IDE-SECIR model:
 *
 * ExponentialDecay expdecay(1.0);
 * StateAgeFunctionWrapper prob(expdecay);
 * model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
 *
 */
struct StateAgeFunctionWrapper {

    /**
     * @brief Constructs a new StateAgeFunctionWrapper object
     * 
     * @param[in] init_function specifies the initial function.
     */
    StateAgeFunctionWrapper(StateAgeFunction& init_function)
        : m_function(init_function.clone())
    {
    }

    /**
     * @brief Copy constructor. 
     */
    StateAgeFunctionWrapper(StateAgeFunctionWrapper const& other)
        : m_function(other.m_function->clone())
    {
    }

    /**
     * @brief Move constructor. 
     */
    StateAgeFunctionWrapper(StateAgeFunctionWrapper&& other) = default;

    /**
     * @brief Copy assignment. 
     */
    StateAgeFunctionWrapper& operator=(StateAgeFunctionWrapper const& other)
    {
        m_function = other.m_function->clone();
        return *this;
    }

    /**
     * @brief Move assignment. 
     */
    StateAgeFunctionWrapper& operator=(StateAgeFunctionWrapper&& other) = default;

    /**
     * @brief Destructor. 
     */
    ~StateAgeFunctionWrapper() = default;

    /**
     * @brief Comparison operator. 
     */
    bool operator==(const StateAgeFunctionWrapper& other) const
    {
        return (m_function->get_state_age_function_type() == other.get_state_age_function_type() &&
                m_function->get_parameter() == other.get_parameter());
    }

    /**
     * @brief Set the StateAgeFunction object
     *
     * @param[in] new_function function that we want to set member m_function to.
     */
    void set_state_age_function(StateAgeFunction& new_function)
    {
        m_function = new_function.clone();
    }

    /**
     * @brief Get type of StateAgeFunction, i.e. which derived class is used.
     *
     * @param[out] string 
     */
    std::string get_state_age_function_type() const
    {
        return m_function->get_state_age_function_type();
    }

    /**
     * @brief Accesses eval of m_function.
     *
     * @param[in] state_age Time at which the function is evaluated.
     * @return Evaluation of the function at state_age. 
     */
    ScalarType eval(ScalarType state_age) const
    {
        return m_function->eval(state_age);
    }

    /**
     * @brief Get the m_parameter object of m_function.
     * 
     * @return ScalarType 
     */
    ScalarType get_parameter() const
    {
        return m_function->get_parameter();
    }

    /**
     * @brief Set the m_parameter object of m_function. 
     * 
     * @param[in] new_parameter New parameter for StateAgeFunction.
     */
    void set_parameter(ScalarType new_parameter)
    {
        m_function->set_parameter(new_parameter);
    }

    ScalarType get_support_max(ScalarType dt, ScalarType tol = 1e-10) const
    {
        return m_function->get_support_max(dt, tol);
    }

private:
    std::unique_ptr<StateAgeFunction> m_function; ///< Stores StateAgeFunction that is used in Wrapper.
};

} // namespace mio

#endif // STATEAGEFUNCTION_H
