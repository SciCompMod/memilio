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
     * @param f rhs of the ODE
     * @param t0 initial point of independent variable t
     * @param y0 value of y at t0
     * @param dt_init initial integration step size
     * @param core implements the solution method
     */
    template <class F, class Vector>
    OdeIntegrator(F&& f, double t0, Vector&& y0, double dt_init, std::shared_ptr<IntegratorCore> core)
        : m_f(std::forward<F>(f))
        , m_result(t0, y0)
        , m_dt(dt_init)
        , m_next_dt(dt_init)
        , m_core(core)
    {
    }

    /**
     * @brief advance the integrator.
     * @param tmax end point. must be greater than get_t().back()
     */
    Eigen::Ref<Eigen::VectorXd> advance(double tmax);

    TimeSeries<double>& get_result()
    {
        return m_result;
    }

    const TimeSeries<double>& get_result() const
    {
        return m_result;
    }

    /**
     * @brief returns the time step width determined by the IntegratorCore for the next integration step
    */
    double get_dt() const
    {
        return m_next_dt;
    }

    void set_integrator(std::shared_ptr<IntegratorCore> integrator)
    {
        m_core = integrator;
    }

private:
    DerivFunction m_f;
    TimeSeries<double> m_result;
    double m_dt;
    double m_next_dt;
    std::shared_ptr<IntegratorCore> m_core;
};

} // namespace mio

#endif // INTEGRATOR_H
