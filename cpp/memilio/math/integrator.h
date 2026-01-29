/* 
* Copyright (C) 2020-2026 MEmilio
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

#include "memilio/config.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include <memory>
#include <functional>

namespace mio
{

/**
 * @brief Function template to be integrated.
 * @tparam FP A floating point type, e.g. double.
 */
template <typename FP>
using DerivFunction =
    std::function<void(Eigen::Ref<const Eigen::VectorX<FP>> y, FP t, Eigen::Ref<Eigen::VectorX<FP>> dydt)>;

/**
 * @brief Interface class defining the integration step used in a SystemIntegrator.
 * @tparam FP A floating point type, e.g. double.
 * @tparam Integrands One or more function types used for defining the right hand side of a system of equations.
 */
template <typename FP, class... Integrands>
class IntegratorCore
{
public:
    /**
     * @brief Initialize an IntegratorCore.
     * Step size bounds are needed for adaptive integrators, see the step method for more detail.
     * Fixed size steppers ignore those bounds and may use the default constructor for FP.
     * @param dt_min Lower bound to the step size dt, as used in the step method.
     * @param dt_max Upper bound to the step size dt, as used in the step method.
     */
    IntegratorCore(const FP& dt_min, const FP& dt_max)
        : m_dt_min(dt_min)
        , m_dt_max(dt_max)
    {
    }

    virtual ~IntegratorCore() {};

    virtual std::unique_ptr<IntegratorCore<FP, Integrands...>> clone() const = 0;

    /**
     * @brief Make a single integration step.
     *
     * The behaviour of this method changes when the integration scheme has adaptive step sizing. 
     * These changes are noted in the parentheses (...) below.
     * Adaptive integrators use the bounds dt_min and dt_max for dt, accessible through the IntegratorCore member
     * functions get_dt_min() and get_dt_max(), respectively. Fixed step integrators ignore these values.
     * The adaptive step sizing is considered to be successful, if a step of at least size dt_min sufficed tolerances.
     * Tolerances are defined in each implementation, usually using a criterion with absolute and relative tolerances.
     * Even if the step sizing failed, the integrator will make a step of at least size dt_min.
     *
     * @param[in] fs One or more function(s) defining the right hand side of the IVP.
     *     May be called multiple times with different arguments.
     * @param[in] yt The known value of y at time t.
     * @param[in,out] t The current time. It will be increased by dt.
     *     (If adaptive, the increment is instead within [dt_min, dt].)
     * @param[in,out] dt The current step size h=dt. Will not be changed.
     *     (If adaptive, the given dt is used as the maximum step size, and must be within [dt_min, dt_max].
     *      During integration, dt is adjusted in [dt_min, dt_max] to have an optimal size for the next step.)
     * @param[out] ytp1 Set to the approximated value of y at time t + dt.
     *     (If adaptive, this time may be smaller, but it is at least t + dt_min, at most t + dt_max.
     *      Note that the increment on t may be different from the returned value of dt.)
     * @return Always true for nonadaptive methods.
     *     (If adaptive, returns whether the adaptive step sizing was successful.)
     */
    virtual bool step(const Integrands&... fs, Eigen::Ref<const Eigen::VectorX<FP>> yt, FP& t, FP& dt,
                      Eigen::Ref<Eigen::VectorX<FP>> ytp1) const = 0;

    /**
     * @brief Access lower bound to the step size dt.
     * These values will only be used by adaptive steppers. Fixed size steppers ignore them.
     * @return A reference to the minimum possible value of dt.
     * @{
     */
    FP& get_dt_min()
    {
        return m_dt_min;
    }
    const FP& get_dt_min() const
    {
        return m_dt_min;
    }
    /** @} */

    /**
     * @brief Access upper bound to the step size dt.
     * These values will only be used by adaptive steppers. Fixed size steppers ignore them.
     * @return A reference to the maximum possible value of dt.
     * @{
     */
    FP& get_dt_max()
    {
        return m_dt_max;
    }
    const FP& get_dt_max() const
    {
        return m_dt_max;
    }
    /** @} */

private:
    FP m_dt_min, m_dt_max; /// Bounds to step size dt.
};

/// @brief Interface for defining solvers for ODE problems. Also @see IntegratorCore and @see CompartmentalModel.
template <class FP>
using OdeIntegratorCore = IntegratorCore<FP, DerivFunction<FP>>;

/// @brief Interface for defining solvers for SDE problems. Also @see IntegratorCore and @see StochasticModel.
template <class FP>
using SdeIntegratorCore = IntegratorCore<FP, DerivFunction<FP>, DerivFunction<FP>>;

/**
 * @brief Integrate a system of equations over time.
 * @tparam FP A floating point type, e.g. double.
 * @tparam Integrands One or more function types used for defining the right hand side of a system of equations.
 * How multiple Integrands are combined depends on the IntegratorCore implementation.
 */
template <typename FP, class... Integrands>
class SystemIntegrator
{
public:
    using Core = IntegratorCore<FP, Integrands...>;
    /**
     * @brief Create an integrator for a specific IVP.
     * @param[in] core Implements the solver.
     */
    SystemIntegrator(std::unique_ptr<Core>&& core)
        : m_core(std::move(core))
        , m_is_adaptive(false)
    {
    }

    SystemIntegrator(const SystemIntegrator& other)
        : m_core(other.m_core->clone())
        , m_is_adaptive(other.m_is_adaptive)
    {
    }

    SystemIntegrator& operator=(const SystemIntegrator& other)
    {
        if (this != &other) {
            m_core        = other.m_core->clone();
            m_is_adaptive = other.m_is_adaptive;
        }
        return *this;
    }

    /**
     * @brief Advance the integrator.
     * @param[in] fs Integrands passed to the integrator, e.g. a wrapper for `get_derivatives`.
     * @param[in] tmax Time end point. Must be greater than results.get_last_time().
     * @param[in, out] dt Initial integration step size. May be changed by the IntegratorCore.
     * @param[in, out] results List of results. Must contain at least one time point. The last entry is used as
     * initial time and value. A new entry is added for each integration step.
     * @return A reference to the last value in the results time series.
     */
    Eigen::Ref<Eigen::VectorX<FP>> advance(const Integrands&... fs, const FP tmax, FP& dt, TimeSeries<FP>& results)
    {
        // hint at std functions for ADL
        using std::ceil;
        using std::fabs;
        using std::max;
        using std::min;
        const FP t0  = results.get_last_time();
        const FP eps = Limits<FP>::zero_tolerance();
        assert(tmax > t0);
        assert(dt > 0);

        const auto num_steps =
            static_cast<size_t>(ceil((tmax - t0) / dt)); // estimated number of time steps (if equidistant)

        results.reserve(results.get_num_time_points() + num_steps);

        bool step_okay = true;

        FP dt_copy; // used to check whether step sizing is adaptive
        FP dt_restore     = 0.0; // used to restore dt, if dt was decreased to reach tmax
        FP dt_min_restore = m_core->get_dt_min(); // used to restore dt_min, if it was decreased to reach tmax
        FP t              = t0;

        for (size_t i = results.get_num_time_points() - 1; fabs((tmax - t) / (tmax - t0)) > eps; ++i) {
            // We don't make time steps too small as the error estimator of an adaptive integrator
            //may not be able to handle it. this is very conservative and maybe unnecessary,
            //but also unlikely to happen. may need to be reevaluated.

            if (dt > tmax - t) {
                dt_restore = dt;
                dt         = tmax - t;
                // if necessary, also reduce minimal step size such that we do not step past tmax
                m_core->get_dt_min() = min<FP>(tmax - t, m_core->get_dt_min());
                // if dt_min was reduced, the following step will be the last due to dt == dt_min (see step method)
                // dt_min must be restored after this loop
            }

            dt_copy = dt;

            results.add_time_point();
            step_okay &= m_core->step(fs..., results[i], t, dt, results[i + 1]);
            results.get_last_time() = t;

            // if dt has been changed by step, register the current m_core as adaptive.
            m_is_adaptive |= !floating_point_equal<FP>(dt, dt_copy, eps);
        }
        m_core->get_dt_min() = dt_min_restore; // restore dt_min
        // if dt was decreased to reach tmax in the last time iteration,
        // we restore it as it is now probably smaller than required for tolerances
        dt = max<FP>(dt, dt_restore);

        if (m_is_adaptive) {
            if (!step_okay) {
                log_warning("Adaptive step sizing failed. Forcing an integration step of size dt_min.");
            }
            else if (fabs((tmax - t) / (tmax - t0)) > eps) {
                log_warning("Last time step too small. Could not reach tmax exactly.");
            }
            else {
                log_info("Adaptive step sizing successful to tolerances.");
            }
        }

        return results.get_last_value();
    }

    /**
     * @brief Change the IntegratorCore used for integration.
     * @param core The new IntegratorCore.
     */
    void set_integrator_core(std::unique_ptr<Core>&& core)
    {
        m_core.swap(core);
        m_is_adaptive = false;
    }

    /**
     * @brief Access the IntegratorCore used for integration.
     * @return A reference to the IntegratorCore used for integration.
     * @{
     */
    Core& get_integrator_core()
    {
        return *m_core;
    }
    const Core& get_integrator_core() const
    {
        return *m_core;
    }

private:
    std::unique_ptr<Core> m_core;
    bool m_is_adaptive;
};

/**
 * @brief Solver for a system of initial value problems (IVPs) consisting of ordinary differential equations (ODEs).
 * The IVPs are of the form y'(t) = f(y(t), t), y(0) = y0.
 */
template <typename FP>
using OdeIntegrator = SystemIntegrator<FP, DerivFunction<FP>>;

/**
 * @brief Solver for a system of initial value problems (IVPs) consisting of stochastic differential equations (SDEs).
 * The IVPs are of the form dY'(t) = f(Y(t), t)dt + g(Y(t), t)dW(t), y(0) = y0.
 * Each summand corresponds to one DerivFunction. 
 */
template <typename FP>
using SdeIntegrator = SystemIntegrator<FP, DerivFunction<FP>, DerivFunction<FP>>;

} // namespace mio

#endif // INTEGRATOR_H
