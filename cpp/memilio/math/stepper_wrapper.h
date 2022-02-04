/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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

#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/math/integrator.h"

// functions and operators neccessary for a Contolled Stepper to work with Eigen::VectorXd
// these have to be declared *before* the includes

namespace std
{
/// @brief elementwise absolute values of an Eigen::VextorXd
Eigen::VectorXd abs(Eigen::VectorXd x);
} // namespace std

/// @brief elementwise addition of a scalar and an Eigen::VextorXd
Eigen::VectorXd operator+(const double s, const Eigen::VectorXd& v);

/// @brief elementwise division of two Eigen::VextorXd
Eigen::VectorXd operator/(const Eigen::VectorXd& v, const Eigen::VectorXd& w);

GCC_CLANG_DIAGNOSTIC(push)
GCC_CLANG_DIAGNOSTIC(ignored "-Wshadow")
GCC_CLANG_DIAGNOSTIC(ignored "-Wlanguage-extension-token")
MSVC_WARNING_DISABLE_PUSH(4127)
#include "boost/numeric/odeint/algebra/vector_space_algebra.hpp"
#include "boost/numeric/odeint/stepper/controlled_runge_kutta.hpp"
#include "boost/numeric/odeint/stepper/runge_kutta4.hpp"
#include "boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp"
MSVC_WARNING_POP
GCC_CLANG_DIAGNOSTIC(pop)

namespace boost
{
namespace numeric
{
    namespace odeint
    {
        // create struct specialization for Eigen::VectorXd of
        // the l-infinity norm used by controlled_runge_kutter
        template <>
        struct vector_space_norm_inf<Eigen::VectorXd> {
            typedef double result_type; // = Eigen::VectorXd::Scalar
            double operator()(Eigen::VectorXd x) const
            {
                return x.lpNorm<Eigen::Infinity>();
            }
        };
    } // namespace odeint
} // namespace numeric
} // namespace boost

// Wrapper implementing IntegratorCore, can be used for several boost::numeric::odeint Controlled Steppers

namespace mio
{

/**
 * @brief Creates and manages an instance of a boost::numeric::odeint::controlled_runge_kutta
 * integrator, wrapped as mio::IntegratorCore.
 */
template <template <class State = Eigen::VectorXd, class Value = double, class Deriv = State, class Time = double,
                    class Algebra    = boost::numeric::odeint::vector_space_algebra,
                    class Operations = typename boost::numeric::odeint::operations_dispatcher<State>::operations_type,
                    class Resizer    = boost::numeric::odeint::never_resizer>
          class ControlledStepper>
class ControlledStepperWrapper : public mio::IntegratorCore
{
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
        : m_abs_tol(abs_tol)
        , m_rel_tol(rel_tol)
        , m_dt_min(dt_min)
        , m_dt_max(dt_max)
        , m_stepper(boost::numeric::odeint::default_error_checker<typename ControlledStepper<>::value_type,
                                                                  typename ControlledStepper<>::algebra_type,
                                                                  typename ControlledStepper<>::operations_type>(
                        abs_tol, rel_tol),
                    boost::numeric::odeint::default_step_adjuster<typename ControlledStepper<>::value_type,
                                                                  typename ControlledStepper<>::time_type>(dt_max))
    // for more options see: boost/boost/numeric/odeint/stepper/controlled_runge_kutta.hpp
    {
    }

    /**
    * @brief Make a single integration step of a system of ODEs and adapt step width
    * @param[in] yt value of y at t, y(t)
    * @param[in,out] t current time step h=dt
    * @param[in,out] dt current time step h=dt
    * @param[out] ytp1 approximated value y(t+1)
    */
    bool step(const mio::DerivFunction& f, Eigen::Ref<Eigen::VectorXd const> yt, double& t, double& dt,
              Eigen::Ref<Eigen::VectorXd> ytp1) const override
    {
        // copy y(t) to dxdt, since we use the scheme try_step(sys, inout, t, dt) with sys=f, inout=y(t) for
        // in-place computation. This is similiar to do_step, but it updates t and dt
        Eigen::VectorXd dxdt = yt;
        const double t_old   = t; // t is updated by try_step on a successfull step
        do {
            m_stepper.try_step(
                // reorder arguments of the DerivFunction f for the stepper
                [&](const Eigen::VectorXd& x, Eigen::VectorXd& dxds, double s) {
                    dxds.resizeLike(x); // try_step calls sys with a vector of size 0 for some reason
                    f(x, s, dxds);
                },
                dxdt, t, dt);
        } while (t == t_old && dt > m_dt_min);
        ytp1 = dxdt;
        return dt > m_dt_min;
    }

    /// @param tol the required absolute tolerance for comparison of the iterative approximation
    void set_abs_tolerance(double abs_tol)
    {
        m_abs_tol = abs_tol;
        set_stepper(); // TODO : replace this call using a custom version of default_error_checker and default_step_adjuster
    }

    /// @param tol the required relative tolerance for comparison of the iterative approximation
    void set_rel_tolerance(double rel_tol)
    {
        m_rel_tol = rel_tol;
        set_stepper(); // TODO : replace this call using a custom version of default_error_checker and default_step_adjuster
    }

    /// @param dt_min sets the minimum step size
    void set_dt_min(double dt_min)
    {
        m_dt_min = dt_min;
    }

    /// @param dt_max sets the maximum step size
    void set_dt_max(double dt_max)
    {
        m_dt_max = dt_max;
        set_stepper(); // TODO : replace this call using a custom version of default_error_checker and default_step_adjuster
    }

private:
    void set_stepper()
    {
        m_stepper = boost::numeric::odeint::controlled_runge_kutta<ControlledStepper<>>(
            boost::numeric::odeint::default_error_checker<typename ControlledStepper<>::value_type,
                                                          typename ControlledStepper<>::algebra_type,
                                                          typename ControlledStepper<>::operations_type>(m_abs_tol,
                                                                                                         m_rel_tol),
            boost::numeric::odeint::default_step_adjuster<typename ControlledStepper<>::value_type,
                                                          typename ControlledStepper<>::time_type>(m_dt_max));
    }
    double m_abs_tol, m_rel_tol, m_dt_min, m_dt_max;
    mutable boost::numeric::odeint::controlled_runge_kutta<ControlledStepper<>> m_stepper;
};

} // namespace mio

#endif