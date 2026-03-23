/* 
* Copyright (C) 2020-2026 MEmilio
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
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/math/smoother.h"
#include "memilio/math/floating_point.h"

#include "boost/math/distributions/gamma.hpp"
#include "boost/math/distributions/lognormal.hpp"
#include "memilio/utils/logging.h"

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
 *     represents a cumulative distribution function. Such functions are also called survival functions.
 *  b) Arbitrary non-negative functions used for parameters such as TransmissionProbabilityOnContact.
 * 
 * Derived classes must implement the eval method which implements the actual function that is evaluated at some state age.
 * This function can depend on the parameter 'scale' to scale the function and on the parameter 'location' to shift the
 * function. Location should be a positive number to fulfill the characteristics of a TransitionDistribution and 
 * scale has to be positive.
 * For a Function F we normally use these parameters at state age x as F(x,location,scale)=F((x-location)/scale). 
 * These two parameters are optional and a derived class does not have to use them.
 * Additionally there is one parameter which specifies the distribution.
 *
 * The derived classes must also implement the clone_impl method which allows to deepcopy the derived class.
 * 
 * The get_support_max method is virtual and implements a basic version to determine the maximum of the support. 
 * For some derived classes there is a more efficient way (see e.g., SmootherCosine) to do this which is 
 * why it can be overridden. The base class implementation uses the fact that the StateAgeFunction is monotonously 
 * decreasing. This is no limitation as the support is only needed for StateAgeFunctions of Type a) as given above.
 * For classes of type b) a dummy implementation logging an error and returning -2 for get_support_max() should be 
 * implemented.
 *
 * The get_mean method is virtual and implements a basic version to determine the mean value of the StateAgeFunction. 
 * The base class implementation uses the fact that the StateAgeFunction is a survival function 
 * (i.e. 1-CDF for any cumulative distribution function CDF). 
 * Therefore, the base class implementation should only be used for StateAgeFunction%s of type a).
 * For some derived classes there is a more efficient way (see e.g., ExponentialSurvivalFunction) to do this which is 
 * why it can be overridden. 
 *
 * See ExponentialSurvivalFunction, SmootherCosine and ConstantFunction for examples of derived classes.
 */
template <typename FP>
struct StateAgeFunction {

    /**
     * @brief Constructs a new StateAgeFunction object
     * 
     * @param[in] init_distribution_parameter Specifies the initial distribution parameter of the function.
     * @param[in] init_location A parameter to shift the function. 
     * @param[in] init_scale A parameter to scale the function. Parameter has to be positive.
     */
    StateAgeFunction(FP init_distribution_parameter, FP init_location = 0, FP init_scale = 1)
        : m_distribution_parameter{init_distribution_parameter}
        , m_location{init_location}
        , m_scale{init_scale}
        , m_mean{-1.0} // Initialize mean as not set.
        , m_support_max{-1.0} // Initialize support maximum as not set.
    {
        if (m_scale <= 0) {
            log_error("The scale parameter of a StateAgeFunction has to be positive. Set scale to 1.");
            m_scale = 1.0;
        }
    }

    /**
     * @brief Virtual destructor.
     */
    virtual ~StateAgeFunction() = default;

    /**
     * @brief Copy constructor.
     */
    StateAgeFunction(const StateAgeFunction<FP>& other) = default;

    /**
     * @brief Move constructor.
     */
    StateAgeFunction(StateAgeFunction<FP>&& other) = default;

    /**
     * @brief Copy assignment operator.
     */
    StateAgeFunction<FP>& operator=(const StateAgeFunction<FP>& other) = default;

    /**
     * @brief Move assignment operator.
     */
    StateAgeFunction<FP>& operator=(StateAgeFunction<FP>&& other) = default;

    /**
     * @brief Comparison operator.
     */
    bool operator==(const StateAgeFunction<FP>& other) const
    {
        return (typeid(*this).name() == typeid(other).name() &&
                m_distribution_parameter == other.get_distribution_parameter() && m_location == other.get_location() &&
                m_scale == other.get_scale());
    }

    /**
     * @brief Here a pure virtual function is defined that depends on the state_age. 
     *
     * The defined function ususally depends on some function parameter.
     * 
     * @param[in] state_age Time at which the function is evaluated.
     */
    virtual FP eval(FP state_age) = 0;

    /**
     * @brief Get the m_distribution_parameter object.
     * 
     * Can be used to access the m_distribution_parameter object, which specifies the used function.
     * 
     * @return ScalarType 
     */
    FP get_distribution_parameter() const
    {
        return m_distribution_parameter;
    }

    /**
     * @brief Set the m_distribution_parameter object.
     * 
     * Can be used to set the m_distribution_parameter object, which specifies the used function.
     * The maximum support of a function may be costly to evaluate. In order to not always reevaluate or recompute the
     * support when the user asks for it, a cached value is used. If m_support_max is set to -1, the cached value is
     * deleted and a recomputation is done the next time the user asks for the support. As the support (potentially)
     * depends on the m_distribution_parameter object, the cached value has to be deleted. For details see get_support_max().
     * The same applies to the m_mean object. See get_mean().
     *
     *@param[in] new_distribution_parameter New parameter for StateAgeFunction.
     */
    void set_distribution_parameter(FP new_distribution_parameter)
    {
        m_distribution_parameter = new_distribution_parameter;
        m_support_max            = -1.;
        m_mean                   = -1.;
    }

    /**
     * @brief Get the m_location object.
     * 
     * Can be used to access the m_location object, which specifies the shift of the function.
     * 
     * @return ScalarType 
     */
    FP get_location() const
    {
        return m_location;
    }

    /**
     * @brief Set the m_location object.
     * 
     * Can be used to set the m_location object, which specifies the shift of the function.
     * The maximum support of a function may be costly to evaluate. In order to not always reevaluate or recompute the
     * support when the user asks for it, a cached value is used. If m_support_max is set to -1, the cached value is
     * deleted and a recomputation is done the next time the user asks for the support. As the support (potentially)
     * depends on the m_location object, the cached value has to be deleted. For details see get_support_max().
     * The same applies to the m_mean object. See get_mean().
     *
     *@param[in] new_location New location for StateAgeFunction.
     */
    void set_location(FP new_location)
    {
        m_location    = new_location;
        m_support_max = -1.;
        m_mean        = -1.;
    }

    /**
     * @brief Get the m_scale object.
     * 
     * Can be used to access the m_scale object, which is used to scale the function.
     * 
     * @return ScalarType 
     */
    FP get_scale() const
    {
        return m_scale;
    }

    /**
     * @brief Set the m_scale object.
     * 
     * Can be used to access the m_scale object, which is used to scale the function.
     * The maximum support of a function may be costly to evaluate. In order to not always reevaluate or recompute the
     * support when the user asks for it, a cached value is used. If m_support_max is set to -1, the cached value is
     * deleted and a recomputation is done the next time the user asks for the support. As the support (potentially)
     * depends on the m_scale object, the cached value has to be deleted. For details see get_support_max().
     * The same applies to the m_mean object. See get_mean().
     *
     *@param[in] new_scale New Scale for StateAgeFunction.
     */
    void set_scale(FP new_scale)
    {
        if (new_scale <= 0) {
            log_error("The scale parameter of a StateAgeFunction has to be positive. Set scale to 1.");
            new_scale = 1;
        }
        m_scale       = new_scale;
        m_support_max = -1.;
        m_mean        = -1.;
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
    virtual FP get_support_max(FP dt, FP tol = 1e-10)
    {
        FP support_max = 0.0;

        if (!floating_point_equal<FP>(m_support_tol, tol, 1e-14) ||
            floating_point_equal<FP>(m_support_max, -1.0, 1e-14)) {
            while (eval(support_max) >= tol) {
                support_max += dt;
            }

            m_support_max = support_max;
            m_support_tol = tol;
        }

        return m_support_max;
    }

    /**
     * @brief Computes the mean value of the function using the time step size dt and some tolerance tol. 
     * 
     * This is a basic version to determine the mean value of a survival function
     * through numerical integration of the integral that describes the expected value.
     * This basic implementation is only valid if the StateAgeFunction is of type a) since we assume that the considered
     * function converges to zero. Otherwise it should be overridden.
     *
     * For some specific derivations of StateAgeFunction%s there are more efficient ways to determine the 
     * mean value which is why this member function is virtual and can be overridden (see, e.g., ExponentialSurvivalFunction).
     * The mean value is only needed for StateAgeFunction%s that are used as TransitionDistribution%s. 
     *
     * @param[in] dt Time step size used for the numerical integration. 
     * @param[in] tol The maximum support used for numerical integration is calculated using this tolerance. 
     * @return ScalarType mean value.
     */
    virtual FP get_mean(FP dt = 1.0, FP tol = 1e-10)
    {
        using std::ceil;
        if (!floating_point_equal<FP>(m_mean_tol, tol, 1e-14) || floating_point_equal<FP>(m_mean, -1., 1e-14)) {
            // Integration using trapezoidal rule.
            // Compute value for i=0.
            FP mean         = 0.5 * dt * eval(FP(0 * dt));
            FP supp_max_idx = ceil(get_support_max(dt, tol) / dt);
            // We start with i=1 since the value for i=0 was already considered above when defining the variable mean.
            // Note that we do not consider indices i>=supp_max_idx since it holds eval(i*dt)<tol for all i>=supp_max_idx
            // by definition of the support_max.
            for (int i = 1; i < supp_max_idx; i++) {
                mean += dt * eval(FP(i * dt));
            }

            m_mean     = mean;
            m_mean_tol = tol;
        }
        return m_mean;
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
    std::unique_ptr<StateAgeFunction<FP>> clone() const
    {
        return std::unique_ptr<StateAgeFunction<FP>>(clone_impl());
    }

protected:
    /**
     * @brief Pure virtual method that implements cloning.
     */
    virtual StateAgeFunction<FP>* clone_impl() const = 0;

    FP m_distribution_parameter; ///< Parameter for function in derived class.
    FP m_location; ///< Location parameter for function in derived class.
    FP m_scale; ///< Scale parameter for function in derived class.
    FP m_mean; ///< Mean value of the function.
    FP m_mean_tol{-1.0}; ///< Tolerance for computation of the mean (initialize as not set).
    FP m_support_max; ///< Maximum of the support of the function.
    FP m_support_tol{-1.0}; ///< Tolerance for computation of the support (initialize as not set).
};

/**************************************
* Derived classes of StateAgeFunction *
***************************************/

/**
 * @brief Class that defines the survival function corresponding to the exponential distribution depending on the state age.
 */

template <typename FP>
struct ExponentialSurvivalFunction : public StateAgeFunction<FP> {

    /**
     * @brief Constructs a new ExponentialSurvivalFunction object.
     * 
     * @param[in] init_distribution_parameter Specifies the initial function parameter of the function.
     * @param[in] init_location Location parameter to shift the ExponentialSurvivalFunction function. 
     *      Should be a positive number to fulfill characteristics of a TransitionDistribution.
     */
    ExponentialSurvivalFunction(FP init_distribution_parameter, FP init_location = 0)
        : StateAgeFunction<FP>(init_distribution_parameter, init_location)
    {
    }

    /**
     * @brief Defines exponential decay function depending on state_age.
     *
     * m_distribution_parameter defines how fast the exponential function decays.
     * 
     * @param[in] state_age Time at which the function is evaluated.
     * @return Evaluation of the function at state_age. 
     */
    FP eval(FP state_age) override
    {
        using std::exp;
        if (state_age <= this->m_location) {
            return 1;
        }
        return exp(-this->m_distribution_parameter * (state_age - this->m_location));
    }

    /**
     * @brief Computes the mean value of the function. 
     * 
     * For the exponential distribution, the mean value is the reciprocal of the distribution parameter and shifted with respect to the location parameter.
     *
     * @param[in] dt Time step size used for the numerical integration (unused for ExponentialSurvivalFunction). 
     * @param[in] tol The maximum support used for numerical integration is calculated using this tolerance 
     *  (unused for ExponentialSurvivalFunction). 
     * @return ScalarType mean value.
     */
    FP get_mean(FP dt = 1.0, FP tol = 1e-10) override
    {
        unused(dt);
        unused(tol);
        return 1.0 / this->m_distribution_parameter + this->m_location;
    }

protected:
    /**
     * @brief Implements clone for ExponentialSurvivalFunction.
     * 
     * @return Pointer to StateAgeFunction.
     */
    StateAgeFunction<FP>* clone_impl() const override
    {
        return new ExponentialSurvivalFunction<FP>(*this);
    }
};

/**
 * @brief Class that defines an smoother_cosine function depending on the state age.
 * This function is a StateAgeFunction of type a) and can therefore be seen as a survival function.
 */
template <typename FP>
struct SmootherCosine : public StateAgeFunction<FP> {

    /**
     * @brief Constructs a new SmootherCosine object.
     *
     * This function is a StateAgeFunction of type a) and can therefore be seen as a survival function.
     *
     * @param[in] init_distribution_parameter specifies the initial parameter of the function.
     * @param[in] init_location Location paramter to shift the SmootherCosine function.
     *      Should be a positive number to fulfill characteristics of a TransitionDistribution.
     */
    SmootherCosine(FP init_distribution_parameter, FP init_location = 0)
        : StateAgeFunction<FP>(init_distribution_parameter, init_location)
    {
    }

    /**
     * @brief Defines smoother cosine function depending on state_age.
     *
     * Used function goes through points (0+m_location,1) and (m_distribution_parameter+m_location,0) and is
     *  interpolated in between using a smoothed cosine function.
     *
     * @param[in] state_age Time at which the function is evaluated.
     * @return Evaluation of the function at state_age.
     */
    FP eval(FP state_age) override
    {
        if (state_age <= this->m_location) {
            return 1.0;
        }
        return smoother_cosine<FP>(state_age - this->m_location, 0.0, this->m_distribution_parameter, 1.0, 0.0);
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
    FP get_support_max(FP dt, FP tol = 1e-10) override
    {
        unused(dt);
        unused(tol);
        this->m_support_max = this->m_distribution_parameter + this->m_location;
        return this->m_support_max;
    }

    /**
     * @brief Computes the mean value of the function.
     *
     * For the associated distribution to SmootherCosine, the mean value is 0.5 * m_distribution_parameter + m_location.
     *
     * @param[in] dt Time step size used for the numerical integration (unused for SmootherCosine).
     * @param[in] tol The maximum support used for numerical integration is calculated using this tolerance
     *  (unused for SmootherCosine).
     * @return ScalarType mean value.
     */
    FP get_mean(FP dt = 1.0, FP tol = 1e-10) override
    {
        unused(dt);
        unused(tol);
        return 0.5 * this->m_distribution_parameter + this->m_location;
    }

protected:
    /**
     * @brief Clones unique pointer to a StateAgeFunction.
     *
     * @return std::unique_ptr<StateAgeFunction> unique pointer to a StateAgeFunction.
     */
    StateAgeFunction<FP>* clone_impl() const override
    {
        return new SmootherCosine<FP>(*this);
    }
};

/**
 * @brief Class that defines an GammaSurvivalFunction function depending on the state age.
 * A survival function is defined as 1 - cumulative density function.
 * GammaSurvivalFunction is derived from StateAgeFunction.
 * The shape parameter of the Gamma function is the parameter of the StateAgeFunction.
 * If shape is an unsigned integer, the Gamma distribution simplifies to an Erlang distribution.
 * Does not support automatic differentiation.
 */
struct GammaSurvivalFunction : public StateAgeFunction<ScalarType> {

    /**
     * @brief Constructs a new GammaSurvivalFunction object.
     *
     * @param[in] init_shape Parameter shape of the GammaSurvivalFunction.
     *  For the Erlang distribution, shape has to be a positive integer.
     *  Choosing shape = 1 leads to an exponential function with parameter 1/scale.
     * @param[in] init_location Location paramter to shift the GammaSurvivalFunction.
     *      Should be a positive number to fulfill characteristics of a TransitionDistribution.
     * @param[in] init_scale Parameter shape of the GammaSurvivalFunction.
     *  Corresponds to the inverse of the rate parameter of a Gamma distribution.
     */
    GammaSurvivalFunction(ScalarType init_shape = 1, ScalarType init_location = 0, ScalarType init_scale = 1)
        : StateAgeFunction<ScalarType>(init_shape, init_location, init_scale)
    {
    }

    /**
     * @brief Defines GammaSurvivalFunction depending on state_age.
     *
     * @param[in] state_age Time at which the function is evaluated.
     * @return Evaluation of the function at state_age.
     */
    ScalarType eval(ScalarType state_age) override
    {
        if (state_age <= this->m_location) {
            return 1;
        }
        boost::math::gamma_distribution<ScalarType, boost::math::policies::policy<>> gamma(
            this->m_distribution_parameter, this->m_scale);
        return boost::math::cdf(boost::math::complement(gamma, state_age - this->m_location));
    }

    /**
     * @brief Computes the mean value of the function.
     *
     * For the gamma distribution, the mean value is m_distribution_parameter*m_scale+m_location,
     * where m_distribution_parameter is the shape parameter.
     *
     * @param[in] dt Time step size used for the numerical integration (unused for GammaSurvivalFunction).
     * @param[in] tol The maximum support used for numerical integration is calculated using this tolerance
     *  (unused for GammaSurvivalFunction).
     * @return ScalarType mean value.
     */
    ScalarType get_mean(ScalarType dt = 1.0, ScalarType tol = 1e-10) override
    {
        unused(dt);
        unused(tol);
        return this->m_distribution_parameter * this->m_scale + this->m_location;
    }

protected:
    /**
     * @brief Implements clone for GammaSurvivalFunction.
     *
     * @return Pointer to StateAgeFunction.
     */
    StateAgeFunction<ScalarType>* clone_impl() const override
    {
        return new GammaSurvivalFunction(*this);
    }
};

/**
 * @brief Class that defines an LognormSurvivalFunction function depending on the state age.
 * A survival function is defined as 1 - cumulative density function.
 * Does not support automatic differentiation.
 */
struct LognormSurvivalFunction : public StateAgeFunction<ScalarType> {

    /**
     * @brief Constructs a new LognormSurvivalFunction object.
     *
     * Location and scale parameters are according to these parameters in the python package scipy.
     *
     * @param[in] init_distribution_parameter Specifies the initial function parameter of the function.
     * @param[in] init_location Location parameter of LognormSurvivalFunction. The parameter can be
     *       used to shift the function. Should be non-negative to fulfill the conditions of a
     *       StateAgeFunction.
     * @param[in] init_scale Scale parameter of LognormSurvivalFunction.
     */
    LognormSurvivalFunction(ScalarType init_distribution_parameter, ScalarType init_location = 0,
                            ScalarType init_scale = 1)
        : StateAgeFunction<ScalarType>(init_distribution_parameter, init_location, init_scale)
    {
    }

    /**
     * @brief Defines the value of the LognormSurvivalFunction depending on state_age.
     *
     * @param[in] state_age Time at which the function is evaluated.
     * @return Evaluation of the function at state_age.
     */
    ScalarType eval(ScalarType state_age) override
    {
        if (state_age < this->m_location) {
            return 1;
        }
        boost::math::lognormal_distribution<ScalarType, boost::math::policies::policy<>> logn(
            0.0, this->m_distribution_parameter);
        return boost::math::cdf(boost::math::complement(logn, (state_age - this->m_location) / this->m_scale));
    }

    // Closed form for the mean value is kind of complex, use default implementation.
    // For testing purposes, a class must exist anyway that uses the default implementation.

protected:
    /**
     * @brief Implements clone for LognormSurvivalFunction.
     *
     * @return Pointer to StateAgeFunction.
     */
    StateAgeFunction<ScalarType>* clone_impl() const override
    {
        return new LognormSurvivalFunction(*this);
    }
};

/**
 * @brief Class that defines a constant function.
 */
template <typename FP>
struct ConstantFunction : public StateAgeFunction<FP> {
    /**
     * @brief Constructs a new ConstantFunction object.
     *
     * @param init_distribution_parameter specifies value of the constant function.
     */
    ConstantFunction(FP init_distribution_parameter)
        : StateAgeFunction<FP>(init_distribution_parameter)
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
    FP eval(FP state_age) override
    {
        unused(state_age);
        return this->m_distribution_parameter;
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
    FP get_support_max(FP dt, FP tol = 1e-10) override
    {
        // In case of a ConstantFunction we would have support_max = infinity
        // This type of function is not suited to be a TransitionDistribution
        // Log error and return -2.

        unused(dt);
        unused(tol);
        this->m_support_max = -2.0;

        log_error("This function is not suited to be a TransitionDistribution. Do not call in case of "
                  "StateAgeFunctions of type b); see documentation of StateAgeFunction Base class.");

        return this->m_support_max;
    }

    /**
     * @brief Computes the mean value of the function.
     *
     * For ConstantFunction, the mean value is the function parameter.
     *
     * @param[in] dt Time step size used for the numerical integration (unused for ConstantFunction).
     * @param[in] tol The maximum support used for numerical integration is calculated using this tolerance (unused for ConstantFunction).
     * @return ScalarType mean value.
     */
    FP get_mean(FP dt = 1.0, FP tol = 1e-10) override
    {
        unused(dt);
        unused(tol);
        log_warning("Attention: This function is not suited to be a TransitionDistribution. Do not call in case of "
                    "StateAgeFunctions of type b); see documentation of StateAgeFunction Base class.");
        return this->m_distribution_parameter;
    }

protected:
    /**
     * @brief Clones unique pointer to a StateAgeFunction.
     *
     * @return std::unique_ptr<StateAgeFunction> unique pointer to a StateAgeFunction
     */
    StateAgeFunction<FP>* clone_impl() const override
    {
        return new ConstantFunction<FP>(*this);
    }
};

/**
 * @brief Class that defines the probability density function corresponding to the Erlang distribution with the parameters shape and scale depending on the state age.
 * Class is needed for the initialization of the subcompartments for LCT model.
 * ErlangDensity is derived from StateAgeFunction.
 * The shape parameter of the Erlang function is the distribution parameter of the StateAgeFunction.
 * Attention: The density is not a survival function and does not have the characteristics of a TransitionDistribution!!
 * The function is of the StateAgeFunction-Type b).
 */
struct ErlangDensity : public StateAgeFunction<ScalarType> {

    /**
     * @brief Constructs a new ErlangDensity object.
     *
     * @param[in] init_shape Parameter shape of the ErlangDensity. For the Erlang distribution, shape has to be a positive integer.
      * @param[in] init_scale Parameter scale of the ErlangDensity. Corresponds to the inverse rate parameter.
     */
    ErlangDensity(unsigned int init_shape, ScalarType init_scale)
        : StateAgeFunction<ScalarType>(init_shape, 0.0, init_scale)
    {
    }

    /**
     * @brief Defines ErlangDensity depending on state_age.
     *
     * Parameters scale and shape are used.
     *
     * @param[in] state_age Time at which the function is evaluated.
     * @return Evaluation of the function at state_age.
     */
    ScalarType eval(ScalarType state_age) override
    {
        using std::exp;
        using std::pow;
        if (state_age < 0) {
            return 0;
        }
        int shape = static_cast<int>(this->m_distribution_parameter);
        return pow(state_age / this->m_scale, shape - 1) /
               (this->m_scale * boost::math::factorial<ScalarType>(shape - 1)) * exp(-state_age / this->m_scale);
    }

    /**
     * @brief Computes the maximum of the support of the function.
     *
     * Calculates the smallest time value t where function(tau)=0 for all tau>t.
     *
     * @param[in] dt Time step size.
     * @param[in] tol Tolerance used for cutting the support if the function value falls below.
     * @return ScalarType support_max
     */
    ScalarType get_support_max(ScalarType dt, ScalarType tol = 1e-10) override
    {
        // We are looking for the smallest time value t where function(tau)=0 for all tau>t. Thus support max is bigger
        // than the mean.
        // We use, that the density is monotonically decreasing for tau>mean here.
        ScalarType mean        = this->m_distribution_parameter * this->m_scale;
        ScalarType support_max = dt * static_cast<int>(mean / dt);

        if (!floating_point_equal<ScalarType>(this->m_support_tol, tol, 1e-14) ||
            floating_point_equal<ScalarType>(this->m_support_max, -1.0, 1e-14)) {
            while (eval(support_max) >= tol) {
                support_max += dt;
            }

            this->m_support_max = support_max;
            this->m_support_tol = tol;
        }

        return this->m_support_max;
    }

    /**
     * @brief Computes the mean value of the function.
     *
     * For the Erlang distribution, the mean value is the m_distribution_parameter * m_scale.
     *
     * @param[in] dt Time step size used for the numerical integration (unused for ErlangDensity).
     * @param[in] tol The maximum support used for numerical integration is calculated using this tolerance (unused for ErlangDensity).
     * @return ScalarType mean value.
     */
    ScalarType get_mean(ScalarType dt = 1.0, ScalarType tol = 1e-10) override
    {
        unused(dt);
        unused(tol);
        log_warning("Attention: This function is not suited to be a TransitionDistribution. Do not call in case of "
                    "StateAgeFunctions of type b); see documentation of StateAgeFunction Base class.");
        return this->m_distribution_parameter * this->m_scale;
    }

protected:
    /**
     * @brief Implements clone for ErlangDensity.
     *
     * @return Pointer to StateAgeFunction.
     */
    StateAgeFunction<ScalarType>* clone_impl() const override
    {
        return new ErlangDensity(*this);
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
 * ExponentialSurvivalFunction exponential(1.0);
 * StateAgeFunctionWrapper prob(exponential);
 * model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
 *
 */
template <typename FP>
struct StateAgeFunctionWrapper {

    StateAgeFunctionWrapper()
        : m_function(mio::SmootherCosine<FP>(1.0).clone())
    {
    }
    /**
     * @brief Constructs a new StateAgeFunctionWrapper object
     *
     * @param[in] init_function specifies the initial function.
     */
    StateAgeFunctionWrapper(StateAgeFunction<FP>& init_function)
        : m_function(init_function.clone())
    {
    }

    /**
     * @brief Copy constructor.
     */
    StateAgeFunctionWrapper(StateAgeFunctionWrapper<FP> const& other)
        : m_function(other.m_function->clone())
    {
    }

    /**
     * @brief Move constructor.
     */
    StateAgeFunctionWrapper(StateAgeFunctionWrapper<FP>&& other) = default;

    /**
     * @brief Copy assignment.
     */
    StateAgeFunctionWrapper<FP>& operator=(StateAgeFunctionWrapper<FP> const& other)
    {
        m_function = other.m_function->clone();
        return *this;
    }

    /**
     * @brief Move assignment.
     */
    StateAgeFunctionWrapper<FP>& operator=(StateAgeFunctionWrapper<FP>&& other) = default;

    /**
     * @brief Destructor.
     */
    ~StateAgeFunctionWrapper() = default;

    /**
     * @brief Comparison operator.
     */
    bool operator==(const StateAgeFunctionWrapper<FP>& other) const
    {
        return (m_function->get_state_age_function_type() == other.get_state_age_function_type() &&
                m_function->get_distribution_parameter() == other.get_distribution_parameter() &&
                m_function->get_location() == other.get_location() && m_function->get_scale() == other.get_scale());
    }

    /**
     * @brief Set the StateAgeFunction object
     *
     * @param[in] new_function function that we want to set member m_function to.
     */
    void set_state_age_function(StateAgeFunction<FP>& new_function)
    {
        m_function = new_function.clone();
    }

    /**
     * @brief Get type of StateAgeFunction, i.e. which derived class is used.
     *
     * @return string
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
    FP eval(FP state_age) const
    {
        return m_function->eval(state_age);
    }

    /**
     * @brief Get the m_distribution_parameter object of m_function.
     *
     * @return ScalarType
     */
    FP get_distribution_parameter() const
    {
        return m_function->get_distribution_parameter();
    }

    /**
     * @brief Set the m_distribution_parameter object of m_function.
     *
     * @param[in] new_distribution_parameter New parameter for StateAgeFunction.
     */
    void set_distribution_parameter(FP new_distribution_parameter)
    {
        m_function->set_distribution_parameter(new_distribution_parameter);
    }

    /**
     * @brief Get the m_location object of m_function.
     *
     * @return ScalarType
     */
    FP get_location() const
    {
        return m_function->get_location();
    }

    /**
     * @brief Set the m_location object of m_function.
     *
     * @param[in] new_location New location for StateAgeFunction.
     */
    void set_location(FP new_location)
    {
        m_function->set_location(new_location);
    }
    /**
     * @brief Get the m_scale object of m_function.
     *
     * @return ScalarType
     */
    FP get_scale() const
    {
        return m_function->get_scale();
    }

    /**
     * @brief Set the m_scale object of m_function.
     *
     * @param[in] new_scale New scale for StateAgeFunction.
     */
    void set_scale(FP new_scale)
    {
        m_function->set_scale(new_scale);
    }

    /**
     * @brief Get the m_support_max object of m_function.
     *
     * @param[in] dt Time step size at which function will be evaluated.
     * @param[in] tol Tolerance used for cutting the support if the function value falls below.
     * @return ScalarType m_support_max
     */
    FP get_support_max(FP dt, FP tol = 1e-10) const
    {
        return m_function->get_support_max(dt, tol);
    }

    /**
     * @brief Get the m_mean object of m_function.
     *
     * @param[in] dt Time step size used for the numerical integration.
     * @param[in] tol The maximum support used for numerical integration is calculated using this tolerance.
     * @return ScalarType m_mean
     */
    FP get_mean(FP dt = 1.0, FP tol = 1e-10) const
    {
        return m_function->get_mean(dt, tol);
    }

private:
    std::unique_ptr<StateAgeFunction<FP>> m_function; ///< Stores StateAgeFunction that is used in Wrapper.
};

} // namespace mio

#endif // STATEAGEFUNCTION_H
