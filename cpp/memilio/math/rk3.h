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
#ifndef RK3_H
#define RK3_H

#include "memilio/config.h"
#include "memilio/math/integrator.h"

namespace mio
{

/**
 * @brief Simple explicit RK3 (Kutta's third-order method) integration for ODE y'(t) = f(t,y)
 * 
 * Scheme:
 *   k1 = f(t, y)
 *   k2 = f(t + h/2, y + h/2 * k1)
 *   k3 = f(t + h, y - h * k1 + 2h * k2)
 *   y(t+h) = y(t) + h/6 * (k1 + 4*k2 + k3)
 * 
 * @tparam FP A floating point type, e.g., ScalarType.
 */
template <typename FP = ScalarType>
class RK3IntegratorCore : public IntegratorCore<FP>
{
public:
    RK3IntegratorCore()
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
        Eigen::VectorX<FP> k2(yt.size());
        Eigen::VectorX<FP> k3(yt.size());
        Eigen::VectorX<FP> y_temp(yt.size());

        // Stage 1: k1 = f(t, y)
        f(yt, t, k1);

        // Stage 2: k2 = f(t + h/2, y + h/2 * k1)
        y_temp = yt + (dt / 2) * k1;
        f(y_temp, t + dt / 2, k2);

        // Stage 3: k3 = f(t + h, y - h * k1 + 2h * k2)
        y_temp = yt - dt * k1 + 2 * dt * k2;
        f(y_temp, t + dt, k3);

        // Final: y(t+h) = y(t) + h/6 * (k1 + 4*k2 + k3)
        ytp1 = yt + (dt / 6) * (k1 + 4 * k2 + k3);
        t += dt;
        return true;
    }
};

} // namespace mio

#endif // RK3_H
