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

namespace mio
{
namespace lsecir
{

/**
 * @brief Class that defines an exponential decay function depending on the state age.
 */
struct ErlangDistribution : public StateAgeFunction {

    /**
     * @brief Constructs a new ExponentialDecay object
     * 
     * @param[in] init_parameter Specifies the initial function parameter of the function.
     */
    ErlangDistribution(ScalarType rate, unsigned int shape)
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
        ScalarType result = 0;
        for (int i = 0; i < m_parameter; i++) {
            result += std::pow(m_rate * state_age, i) / (boost::math::factorial<ScalarType>(i));
        }
        return result * std::exp(-m_rate * state_age);
    }

protected:
    /**
     * @brief Implements clone for ExponentialDecay.
     * 
     * @return Pointer to StateAgeFunction.
     */
    StateAgeFunction* clone_impl() const override
    {
        return new ErlangDistribution(*this);
    }

    ScalarType m_rate;
};

Eigen::VectorXd compute_compartments(TimeSeries<ScalarType>&& init, InfectionStateBase Base, InfectionState InfState,
                                     ScalarType gamma, Eigen::Index idx_IncomingFlow, ScalarType dt);

} // namespace lsecir
} // namespace mio

#endif