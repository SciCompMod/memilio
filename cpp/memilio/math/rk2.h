/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
*
* Contact: Martin J. Kuehn@DLR.de
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
#ifndef RK2_H
#define RK2_H

#include "memilio/config.h"
#include "memilio/math/integrator.h"

namespace mio
{

/**
 * @brief Simple explicit RK2 (Midpoint method) integration for ODE y'(t) = f(t,y)
 * 
 * Scheme:
 *   k1 = f(t, y)
 *   k2 = f(t + h/2, y + h/2 * k1)
 *   y(t+h) = y(t) + h * k2
 * 
 * @tparam FP A floating point type, e.g., ScalarType.
 */
template <typename FP = ScalarType>
class RK2IntegratorCore : public IntegratorCore<FP>
{
public:
    RK2IntegratorCore()
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
        Eigen::VectorX<FP> k1(yt.size());
        Eigen::VectorX<FP> y_mid(yt.size());

        // Stage 1: k1 = f(t, y)
        f(yt, t, k1);

        // Stage 2: k2 = f(t + h/2, y + h/2 * k1)
        y_mid = yt + (dt / 2) * k1;
        f(y_mid, t + dt / 2, ytp1); // Store k2 in ytp1

        // Final: y(t+h) = y(t) + h * k2
        ytp1 = yt + dt * ytp1;
        t += dt;
        return true;
    }
};

} // namespace mio

#endif // RK2_H
