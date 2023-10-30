/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "memilio/utils/time_series.h"

#include "memilio/math/eigen.h"
#include <memory>
#include <vector>
#include <functional>
#include <algorithm>

namespace mio
{

/**
 * Function template to be integrated
 */
using DerivFunction =
    std::function<void(Eigen::Ref<const Eigen::VectorXd> y, double t, Eigen::Ref<Eigen::VectorXd> dydt)>;

class IntegratorCore
{
public:
    virtual ~IntegratorCore(){};

    /**
     * @brief Step of the integration with possible adaptive with
     *
     * @param[in] f Right hand side of ODE
     * @param[in] yt value of y at t, y(t)
     * @param[in,out] t current time step h=dt
     * @param[in,out] dt current time step h=dt
     * @param[out] ytp1 approximated value y(t+1)
     */
    virtual bool step(const DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
                      Eigen::Ref<Eigen::VectorXd> ytp1) const = 0;
};

/**
 * Integrate initial value problems (IVP) of ordinary differential equations (ODE) of the form y' = f(y, t), y(t0) = y0.
 */
class OdeIntegrator
{
public:
    /**
     * @brief create an integrator for a specific IVP
     * @param[in] core implements the solution method
     */
    OdeIntegrator(std::shared_ptr<IntegratorCore> core)
        : m_core(core)
    {
    }

    /**
     * @brief Advance the integrator.
     * @param[in] f The rhs of the ODE.
     * @param[in] tmax Time end point. Must be greater than results.get_last_time().
     * @param[in, out] dt Initial integration step size. May be changed by the IntegratorCore.
     * @param[in, out] results List of results. Must contain at least one time point. The last entry is used as
     * intitial time and value. A new entry is added for each integration step.
     * @return The last result.
     */
    Eigen::Ref<Eigen::VectorXd> advance(const DerivFunction& f, const double tmax, double& dt,
                                        TimeSeries<double>& results);

    void set_integrator(std::shared_ptr<IntegratorCore> integrator)
    {
        m_core = integrator;
    }

private:
    std::shared_ptr<IntegratorCore> m_core;
};

} // namespace mio

#endif // INTEGRATOR_H
