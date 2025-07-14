/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Martin J. Kuehn, Daniel Abele
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
#ifndef EULER_H
#define EULER_H

#include "memilio/config.h"
#include "memilio/math/integrator.h"

namespace mio
{

/**
 * @brief Simple explicit euler integration y(t+1) = y(t) + h*f(t,y) for ODE y'(t) = f(t,y)
 * @tparam FP A floating point type, e.g., ScalarType.
 */
template <typename FP>
class EulerIntegratorCore : public IntegratorCore<FP>
{
public:
    EulerIntegratorCore()
        : IntegratorCore<FP>(FP{}, FP{})
    {
    }

    /**
     * @brief Fixed step width of the integration
     *
     * @param[in] yt value of y at t, y(t)
     * @param[in,out] t current time step h=dt
     * @param[in,out] dt current time step h=dt
     * @param[out] ytp1 approximated value y(t+1)
     */
    bool step(const DerivFunction<FP>& f, Eigen::Ref<const Eigen::VectorX<FP>> yt, FP& t, FP& dt,
              Eigen::Ref<Eigen::VectorX<FP>> ytp1) const override
    {
        // we are misusing the next step y as temporary space to store the derivative
        f(yt, t, ytp1);
        ytp1 = yt + dt * ytp1;
        t += dt;
        return true;
    }
};

} // namespace mio

#endif // EULER_H
