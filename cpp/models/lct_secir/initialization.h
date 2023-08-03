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
 * @brief Class that defines an exponential decay function depending on the state age.
 */
struct ErlangDensity : public StateAgeFunction {

    /**
     * @brief Constructs a new ExponentialDecay object
     * 
     * @param[in] init_parameter Specifies the initial function parameter of the function.
     */
    ErlangDensity(ScalarType rate, unsigned int shape)
        : StateAgeFunction(shape)
    {
        m_rate = rate;
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
        if (state_age <= 0) {
            return 0;
        }
        int shape = (int)m_parameter;
        return m_rate * std::pow(m_rate * state_age, shape - 1) / (boost::math::factorial<ScalarType>(shape - 1)) *
               std::exp(-m_rate * state_age);
    }

protected:
    /**
     * @brief Implements clone for ExponentialDecay.
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
 * @brief Class that defines an exponential decay function depending on the state age.
 */
struct ErlangSurvivalFunction : public StateAgeFunction {

    /**
     * @brief Constructs a new ExponentialDecay object
     * 
     * @param[in] init_parameter Specifies the initial function parameter of the function.
     */
    ErlangSurvivalFunction(ScalarType rate, unsigned int shape)
        : StateAgeFunction(shape)
    {
        m_rate = rate;
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
     * @brief Implements clone for ExponentialDecay.
     * 
     * @return Pointer to StateAgeFunction.
     */
    StateAgeFunction* clone_impl() const override
    {
        return new ErlangSurvivalFunction(*this);
    }

    ScalarType m_rate;
};

class Initializer
{
public:
    using ParameterSet = Parameters;

    Initializer(TimeSeries<ScalarType>&& flows, ScalarType dt, InfectionState infectionState_init = InfectionState(),
                ParameterSet&& parameterSet_init = ParameterSet())
        : parameters(parameterSet_init)
        , infectionStates(infectionState_init)
        , m_flows(std::move(flows))
        , m_dt(dt)
    {
    }

    Eigen::VectorXd compute_initializationvector(ScalarType total_population, ScalarType deaths);

    ParameterSet parameters{}; ///< ParameterSet of Model Parameters.
    InfectionState infectionStates; ///< InfectionState specifies number of subcompartments.
private:
    Eigen::VectorXd compute_compartment(InfectionStateBase base, Eigen::Index idx_Incoming_flow,
                                        ScalarType transition_rate);

    ScalarType compute_susceptibles(ScalarType total_population, ScalarType deaths);
    TimeSeries<ScalarType> m_flows;
    ScalarType m_dt{};
};

} // namespace lsecir
} // namespace mio

#endif