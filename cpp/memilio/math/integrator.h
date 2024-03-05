/* 
* Copyright (C) 2020-2024 MEmilio
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

#include <memory>
#include <functional>

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
     * @brief Make a single integration step.
     *
     * The behaviour of this method changes slightly when the integration scheme has adaptive step sizing. 
     * These changes are noted in (brackets) below.
     * Adaptive integrators must have bounds dt_min and dt_max for dt.
     * The adaptive step sizing is considered to have failed, if dt would be adapted to a value strictly below dt_min.
     * Even if the step sizing failed, the integrator must make a step of at least size dt_min.
     *
     * @param[in] f Right hand side of ODE. May be called multiple times with different arguments.
     * @param[in] yt The known value of y at time t.
     * @param[in,out] t The current time. It will be increased by dt.
     *     (If adaptive, the increment is from [dt_min, dt].)
     * @param[in,out] dt The current step size h=dt.
     *     (If adaptive, dt is set to an estimate in [dt_min, dt_max] for the optimal size of the next step.)
     * @param[out] ytp1 Set to the approximated value of y at time t + dt.
     *     (If adaptive, this time may be smaller, but it is at least t + dt_min.
     *      Note that the increment on t may be different from the returned value of dt.)
     * @return Always true.
     *     (If adaptive, returns whether the adaptive step sizing was successfull.)
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
     * @return A reference to the last value in the results time series.
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
