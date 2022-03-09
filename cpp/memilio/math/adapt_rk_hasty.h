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
#ifndef adapt_rk_hasty_H_
#define adapt_rk_hasty_H_

#include "memilio/math/adapt_rk.h"

namespace mio
{

/**
 * @brief Two scheme Runge-Kutta numerical integrator with adaptive step width; faster, but no gurantee for convergence
 *
 * This class is a variation of mio::RKIntegratorCore. Instead of using an elementwise convergence criterion in the
 * step method, only the maximum error is compared against the maximum (mixed) tolerance of all equations. Since
 * these maxima can come from different equations, this is not a valid criterion and cannot guarantee convergence.
 *
 * However, in practice this criterion seems to suffice e.g. for the secir models, and it
 * can be about twice as fast as the regular RKIntegratorCore.
 */
class HastyRKIntegratorCore : public RKIntegratorCore
{
public:
    /// @brief Set up the integrator
    HastyRKIntegratorCore()
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
    HastyRKIntegratorCore(const double abs_tol, const double rel_tol, const double dt_min, const double dt_max)
        : RKIntegratorCore(abs_tol, rel_tol, dt_min, dt_max)
    {
    }
    /**
     * @brief Make a single integration step of a system of ODEs and adapt the step size
     * @param[in] yt value of y at t, y(t)
     * @param[in,out] t current time
     * @param[in,out] dt current time step size h=dt
     * @param[out] ytp1 approximated value y(t+1)
     * 
     * This is variation of mio::RKIntegratorCore::step cannot guarantee convergence but is faster.
     */
    bool step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
              Eigen::Ref<Eigen::VectorXd> ytp1) const override;
};

} // namespace mio

#endif