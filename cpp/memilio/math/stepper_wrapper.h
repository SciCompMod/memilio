/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding
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
#ifndef STEPPER_WRAPPER_H_
#define STEPPER_WRAPPER_H_

#include "memilio/math/integrator.h"
#include "memilio/utils/logging.h"

#include "boost/numeric/odeint/external/eigen/eigen_algebra.hpp"
#include "boost/numeric/odeint/external/eigen/eigen_resize.hpp"
#include "boost/numeric/odeint/stepper/controlled_runge_kutta.hpp"
#include "boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
// #include "boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp" // TODO: reenable once boost bug is fixed

namespace mio
{

/**
 * @brief This is an adaptive IntegratorCore. It creates and manages an instance of a
 * boost::numeric::odeint::controlled_runge_kutta integrator, wrapped as mio::IntegratorCore.
 */
template <typename FP,
          template <class State, class Value, class Deriv, class Time, class Algebra, class Operations, class Resizer>
          class ControlledStepper>
class ControlledStepperWrapper : public mio::IntegratorCore<FP>
{
    using Stepper = boost::numeric::odeint::controlled_runge_kutta<
        ControlledStepper<Eigen::VectorX<FP>, FP, Eigen::VectorX<FP>, FP, boost::numeric::odeint::vector_space_algebra,
                          typename boost::numeric::odeint::operations_dispatcher<Eigen::VectorX<FP>>::operations_type,
                          boost::numeric::odeint::initially_resizer>>;
    static constexpr bool is_fsal_stepper = std::is_same_v<typename Stepper::stepper_type::stepper_category,
                                                           boost::numeric::odeint::explicit_error_stepper_fsal_tag>;
    static_assert(!is_fsal_stepper,
                  "FSAL steppers cannot be used until https://github.com/boostorg/odeint/issues/72 is resolved.");

public:
    /**
     * @brief Set up the integrator
     * @param abs_tol absolute tolerance
     * @param rel_tol relative tolerance 
     * @param dt_min lower bound for time step dt
     * @param dt_max upper bound for time step dt
     */
    ControlledStepperWrapper(double abs_tol = 1e-10, double rel_tol = 1e-5,
                             double dt_min = std::numeric_limits<double>::min(),
                             double dt_max = std::numeric_limits<double>::max())
        : IntegratorCore<FP>(dt_min, dt_max)
        , m_abs_tol(abs_tol)
        , m_rel_tol(rel_tol)
        , m_stepper(create_stepper())
    {
    }

    /**
     * @brief Make a single integration step on a system of ODEs and adapt the step size dt.
     *
     * @param[in] yt Value of y at t, y(t).
     * @param[in,out] t Current time. Will be set to t' in [t+dt_min, t+dt].
     * @param[in,out] dt Current time step size h=dt. Overwritten by an estimated optimal step size for the next step.
     * @param[out] ytp1 The approximated value of y(t').
     */
    bool step(const mio::DerivFunction<FP>& f, Eigen::Ref<Eigen::VectorX<FP> const> yt, FP& t, FP& dt,
              Eigen::Ref<Eigen::VectorX<FP>> ytp1) const override
    {
        using boost::numeric::odeint::fail;
        using std::max;
        assert(0 <= this->get_dt_min());
        assert(this->get_dt_min() <= this->get_dt_max());

        if (dt < this->get_dt_min() || dt > this->get_dt_max()) {
            mio::log_warning("IntegratorCore: Restricting given step size dt = {} to [{}, {}].", dt, this->get_dt_min(),
                             this->get_dt_max());
        }
        // set initial values for exit conditions
        auto step_result = fail;
        bool is_dt_valid = true;
        // copy vectors from the references, since the stepper cannot (trivially) handle Eigen::Ref
        m_ytp1 = ytp1; // y(t')
        m_yt   = yt; // y(t)
        // make a integration step, adapting dt to a possibly larger value on success,
        // or a strictly smaller value on fail.
        // stop only on a successful step or a failed step size adaption (w.r.t. the minimal step size dt_min)
        while (step_result == fail && is_dt_valid) {
            if (dt < this->get_dt_min()) {
                is_dt_valid = false;
                dt          = this->get_dt_min();
            }
            // we use the scheme try_step(sys, in, t, out, dt) with sys=f, in=y(t), out=y(t').
            // this is similiar to do_step, but it can adapt the step size dt. If successful, it also updates t.
            // Note: the resizer used by m_stepper restricts dt to dt_max (via making a failed step)
            if constexpr (!is_fsal_stepper) { // prevent compile time errors with fsal steppers
                step_result = m_stepper.try_step(
                    // reorder arguments of the DerivFunction f for the stepper
                    [&](const Eigen::VectorX<FP>& x, Eigen::VectorX<FP>& dxds, FP s) {
                        f(x, s, dxds);
                    },
                    m_yt, t, m_ytp1, dt);
            }
        }
        // output the new value by copying it back to the reference
        ytp1 = m_ytp1;
        // bound dt from below
        // the last adaptive step (successful or not) may have calculated a new step size smaller than m_dt_min

        dt = max(dt, this->get_dt_min());
        // check whether the last step failed (which means that m_dt_min was still too large to suffice tolerances)
        if (step_result == fail) {
            // adaptive stepping failed, but we still return the result of the last attempt
            t += this->get_dt_min();
            return false;
        }
        else {
            // successfully made an integration step
            return true;
        }
    }

    /// @param tol the required absolute tolerance for comparison of the iterative approximation
    void set_abs_tolerance(FP abs_tol)
    {
        m_abs_tol = abs_tol;
        m_stepper = create_stepper();
    }

    /// @param tol the required relative tolerance for comparison of the iterative approximation
    void set_rel_tolerance(FP rel_tol)
    {
        m_rel_tol = rel_tol;
        m_stepper = create_stepper();
    }

    /// @param dt_min sets the minimum step size
    void set_dt_min(FP dt_min)
    {
        this->get_dt_min() = dt_min;
    }

    /// @param dt_max sets the maximum step size
    void set_dt_max(FP dt_max)
    {
        this->get_dt_max() = dt_max;
        m_stepper          = create_stepper();
    }

private:
    /// @brief (Re)initialize the internal stepper.
    Stepper create_stepper()
    {
        // for more options see: boost/boost/numeric/odeint/stepper/controlled_runge_kutta.hpp
        return Stepper(typename Stepper::error_checker_type(m_abs_tol, m_rel_tol),
                       typename Stepper::step_adjuster_type(this->get_dt_max()));
    }

    FP m_abs_tol, m_rel_tol; ///< Absolute and relative tolerances for integration.
    mutable Eigen::VectorX<FP> m_ytp1, m_yt; ///< Temporary storage to avoid allocations in step function.
    mutable Stepper m_stepper; ///< A stepper instance used for integration.
};

/**
 * @brief This is a non-adaptive IntegratorCore. It creates and manages an instance of an explicit stepper from
 * boost::numeric::odeint, wrapped as mio::IntegratorCore.
 */
template <typename FP,
          template <class State, class Value, class Deriv, class Time, class Algebra, class Operations, class Resizer>
          class ExplicitStepper>
class ExplicitStepperWrapper : public mio::IntegratorCore<FP>
{
public:
    using Stepper =
        ExplicitStepper<Eigen::VectorX<FP>, FP, Eigen::VectorX<FP>, FP, boost::numeric::odeint::vector_space_algebra,
                        typename boost::numeric::odeint::operations_dispatcher<Eigen::VectorX<FP>>::operations_type,
                        boost::numeric::odeint::initially_resizer>;

    /**
     * @brief Set up the integrator.
     */
    ExplicitStepperWrapper()
        : mio::IntegratorCore<FP>(FP{}, FP{})
    {
    }

    /**
     * @brief Make a single integration step on a system of ODEs with fixed step size dt.
     *
     * @param[in] yt Value of y at t, y(t).
     * @param[in,out] t Current time. Overwritten with t+dt.
     * @param[in] dt Current time step size h=dt.
     * @param[out] ytp1 The approximated value of y(t+dt).
     */
    bool step(const mio::DerivFunction<FP>& f, Eigen::Ref<Eigen::VectorX<FP> const> yt, FP& t, FP& dt,
              Eigen::Ref<Eigen::VectorX<FP>> ytp1) const override
    {
        // copy the values from y(t) to ytp1, since we use the scheme do_step(sys, inout, t, dt) with
        // sys=f, inout=y(t) for in-place computation - also, this form is shared by several steppers in boost
        ytp1 = yt;
        m_stepper.do_step(
            // reorder arguments of the DerivFunction f for the stepper
            [&](const Eigen::VectorX<FP>& x, Eigen::VectorX<FP>& dxds, FP s) {
                f(x, s, dxds);
            },
            ytp1, t, dt);
        // update time (it is not modified by do_step)
        t += dt;
        return true; // no step size adaption
    }

private:
    mutable Stepper m_stepper; ///< A stepper instance used for integration.
};

} // namespace mio

#endif
