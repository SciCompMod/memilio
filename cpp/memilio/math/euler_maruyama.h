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
#ifndef MIO_MATH_EULER_MARUYAMA_H
#define MIO_MATH_EULER_MARUYAMA_H

#include "memilio/math/integrator.h"
#include "memilio/math/eigen.h"
#include "memilio/math/utils.h"
#include "memilio/utils/logging.h"

namespace mio
{

/**
 * @brief Simple explicit euler integration y(t+1) = y(t) + h*f(t,y) for ODE y'(t) = f(t,y)
 * @tparam FP A floating point type, e.g., ScalarType.
 */
template <typename FP>
class EulerMaruyamaIntegratorCore : public SdeIntegratorCore<FP>
{
public:
    EulerMaruyamaIntegratorCore()
        : SdeIntegratorCore<FP>(FP{}, FP{})
    {
    }

    /**
     * @brief Fixed step width integrator for stochastic models.
     * Expects population-based integrands, i.e. not flows. 
     * @param[in] f Deterministic part of the stochastic differential equation.
     * @param[in] noise_f Stochastic part of the stochastic differential equation.
     * @param[in] yt Value of y at t, y(t).
     * @param[in,out] t Current time. Overwritten with t+dt.
     * @param[in] dt Current time step size h=dt.
     * @param[out] ytp1 The approximated value of y(t+dt).
     */
    bool step(const DerivFunction<FP>& f, const DerivFunction<FP>& noise_f, Eigen::Ref<const Eigen::VectorX<FP>> yt,
              FP& t, FP& dt, Eigen::Ref<Eigen::VectorX<FP>> ytp1) const override
    {
        using std::sqrt;

        f(yt, t, ytp1);

        Eigen::VectorX<FP> noise = Eigen::VectorX<FP>::Zero(yt.size());
        noise_f(yt, t, noise);

        ytp1 = yt + dt * ytp1 + sqrt(dt) * noise;

        if (!map_to_nonnegative(ytp1)) [[unlikely]] {
            mio::log_error("Failed to rescale values in Euler Maruyama step, total population is negative. "
                           "Redistributing absolute values evenly.");
        }

        t += dt;
        return true;
    }
};

} // namespace mio

#endif // MIO_MATH_EULER_MARUYAMA_H
