/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

#ifndef LCTSECIR_INITIALIZATION_H
#define LCTSECIR_INITIALIZATION_H

#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include <boost/math/special_functions/factorials.hpp>
#include "lct_secir/infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "lct_secir/parameters.h"

namespace mio
{
namespace lsecir
{

/**
 * @brief Class that defines an Erlang density function with the parameters rate and shape depending on the state age.
 * Class is needed for the initialization of the subcompartments for LCT model.
 * ErlangDensity is derived from StateAgeFunction but implements an additional parameter. 
 * The shape parameter of the Erlang function is the parameter of the Stateagefunction and can be changed via the methods of StateAgeFunction. 
 * The rate parameter is not designed to be changed.
 */
struct ErlangDensity : public StateAgeFunction {

    /**
     * @brief Constructs a new ErlangDensity object.
     * 
     * @param[in] rate Parameter rate of the ErlangDensity.
     * @param[in] shape Parameter shape of the ErlangDensity. For the Erlang distribution, shape has to be a positive integer.
     */
    ErlangDensity(ScalarType rate, unsigned int shape)
        : StateAgeFunction(shape)
    {
        m_rate = rate;
    }

    /**
     * @brief Defines ErlangDensity depending on state_age.
     *
     * Parameters rate and shape are used.
     * 
     * @param[in] state_age Time at which the function is evaluated.
     * @return Evaluation of the function at state_age. 
     */
    ScalarType eval(ScalarType state_age) override
    {
        if (state_age <= 0) {
            return 0;
        }
        int shape = (int)m_parameter;
        return m_rate * std::pow(m_rate * state_age, shape - 1) / (boost::math::factorial<ScalarType>(shape - 1)) *
               std::exp(-m_rate * state_age);
    }

protected:
    /**
     * @brief Implements clone for ErlangDensity.
     * 
     * @return Pointer to StateAgeFunction.
     */
    StateAgeFunction* clone_impl() const override
    {
        return new ErlangDensity(*this);
    }

    ScalarType m_rate;
};

/**
 * @brief Class that defines an ErlangSurvivalFunction function depending on the state age.
 * Class is needed for the initialization of the Susceptible compartment from flows for theLCT model.
 * A survival function is defined as 1 - cumulative density function.
 * ErlangSurvivalFunction is derived from StateAgeFunction but implements an additional parameter. 
 * The shape parameter of the Erlang function is the parameter of the Stateagefunction and can be changed via the methods of StateAgeFunction. 
 * The rate parameter is not designed to be changed.
 */
struct ErlangSurvivalFunction : public StateAgeFunction {

    /**
     * @brief Constructs a new ErlangSurvivalFunction object.
     *
     * @param[in] rate Parameter rate of the ErlangSurvivalFunction.
     * @param[in] shape Parameter shape of the ErlangSurvivalFunction. For the Erlang distribution, shape has to be a positive integer.
     */
    ErlangSurvivalFunction(ScalarType rate, unsigned int shape)
        : StateAgeFunction(shape)
    {
        m_rate = rate;
    }

    /**
     * @brief Defines ErlangSurvivalFunction depending on state_age.
     *
     * Parameters rate and shape are used.
     * 
     * @param[in] state_age Time at which the function is evaluated.
     * @return Evaluation of the function at state_age. 
     */
    ScalarType eval(ScalarType state_age) override
    {
        if (state_age <= 0) {
            return 1;
        }
        int result = 0;
        for (int j = 1; j < (int)m_parameter + 1; j++) {
            result += m_rate * std::pow(m_rate * state_age, j - 1) / (boost::math::factorial<ScalarType>(j - 1)) *
                      std::exp(-m_rate * state_age);
        }
        return (1 / m_rate) * result;
    }

protected:
    /**
     * @brief Implements clone for ErlangSurvivalFunction.
     * 
     * @return Pointer to StateAgeFunction.
     */
    StateAgeFunction* clone_impl() const override
    {
        return new ErlangSurvivalFunction(*this);
    }

    ScalarType m_rate;
};

/**
 * @brief Class that can be used to compute an initialization vector out of flows for a LCT Model.
 */
class Initializer
{
public:
    using ParameterSet = Parameters;

    /**
     * @brief Constructs a new Initializer object.
     *
     * @param[in] flows Initalizing TimeSeries with flows fitting to these defined in InfectionTransition. 
     *      Timesteps should be equidistant.
     * @param[in, out] infectionState_init InfectionStates for the Initializer, specifies number of Subcompartments for each infection state.
     * @param[in, out] parameterSet_init Specifies Parameters necessary for the Initializer.
     */
    Initializer(TimeSeries<ScalarType>&& flows, InfectionState infectionState_init = InfectionState(),
                ParameterSet&& parameterSet_init = ParameterSet())
        : parameters(parameterSet_init)
        , infectionStates(infectionState_init)
        , m_flows(std::move(flows))

    {
        m_dt = m_flows.get_time(1) - m_flows.get_time(0);
    }

    /**
     * @brief Core function of Initializer.
     *
     * Computes a vector that can be used for the initalization of an LCT model with the number of persons for each subcompartment.
     *
     * @param[in] total_population The total size of the considered population.
     * @param[in] deaths Number of deceased people from the disease at time 0.
     * @return Vector with a possible initialization for an LCT model computed out of the flows.
     */
    Eigen::VectorXd compute_initializationvector(ScalarType total_population, ScalarType deaths) const;

    ParameterSet parameters{}; ///< ParameterSet with the Initalizers Parameters.
    InfectionState infectionStates; ///< InfectionState specifies number of subcompartments.

private:
    /**
     * @brief Checks constraints of the Initializer inclusive check for parameters.
     */
    void check_constraints() const;

    /**
     * @brief Computes a vector with the number of people in each compartment for one InfectionStateBase.
     *
     * With this function, partial result of compute_initializationvector are achieved.
     *
     * @param[in] base The InfectionStateBase for which the partial result should be calculated.
     * @param[in] idx_incoming_flow Index of the flow which is relevant for the calculation, so the flow to the InfectionStateBase base.
     * @param[in] transition_rate Specifies the transition rate of the InfectionsStateBase. Is equal to 1 / (expected Time in base).
     * @return Vector with a possible initialization for the subcompartments of base.
     */
    Eigen::VectorXd compute_compartment(InfectionStateBase base, Eigen::Index idx_incoming_flow,
                                        ScalarType transition_rate) const;

    /**
     * @brief Computes a vector with the number of people in each compartment for one InfectionStateBase.
     *
     * With this function, a partial result of compute_initializationvector are achieved, namely the value for Susceptibles.
     *
     * @param[in] base The InfectionStateBase for which the partial result should be calculated.
     * @param[in] idx_incoming_flow Index of the flow which is relevant for the calculation, so the flow to the InfectionStateBase base.
     * @param[in] transition_rate Specifies the transition rate of the InfectionsStateBase. Is equal to 1 / (expected Time in base).
     * @return Vector with a possible initialization for the subcompartments of base.
     */
    ScalarType compute_susceptibles(ScalarType total_population, ScalarType deaths) const;

    TimeSeries<ScalarType> m_flows; ///< TimeSeries with the flows which are used to calculate the initial vector.
    ScalarType m_dt{}; ///< Step size of the times in m_flows and time step for the approximation of the integral.
};

} // namespace lsecir
} // namespace mio

#endif