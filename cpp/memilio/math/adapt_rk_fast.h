/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Martin J. Kuehn, Daniel Abele, Rene Schmieding
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
#ifndef ADAPT_RK_FAST_H_
#define ADAPT_RK_FAST_H_

#include "memilio/math/adapt_rk.h"

namespace mio
{

class RKFastIntegratorCore : public RKIntegratorCore
{
public:
    /// @brief Set up the integrator
    RKFastIntegratorCore()
        : RKIntegratorCore()
    {
    }
    /**
     * @brief Set up the integrator
     * @param abs_tol absolute tolerance
     * @param rel_tol relative tolerance 
     * @param dt_min lower bound for time step dt
     * @param dt_max upper bound for time step dt
     */
    RKFastIntegratorCore(const double abs_tol, const double rel_tol, const double dt_min, const double dt_max)
        : RKIntegratorCore(abs_tol, rel_tol, dt_min, dt_max)
    {
    }
    /**
    * @brief Make a single integration step of a system of ODEs and adapt step width
    * @param[in] yt value of y at t, y(t)
    * @param[in,out] t current time step h=dt
    * @param[in,out] dt current time step h=dt
    * @param[out] ytp1 approximated value y(t+1)
    */
    bool step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
              Eigen::Ref<Eigen::VectorXd> ytp1) const override;
};

} // namespace mio

#endif