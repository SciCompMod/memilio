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

#include "memilio/math/integrator.h"

// functions and operators neccessary for a Contolled Stepper to work with Eigen::VectorXd
// these have to be declared *before* the includes

namespace std {
Eigen::VectorXd abs(Eigen::VectorXd x) {
    // elementwise operations are defined on arrays within Eigen
    // casts to and from array supposedly cost no runtime when using compiler optimisation 
    return x.array().abs().matrix();
}
}

#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include <boost/numeric/odeint/stepper/controlled_runge_kutta.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta4.hpp>
#include <boost/numeric/odeint/stepper/adams_bashforth_moulton.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_fehlberg78.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp>
#include <boost/numeric/odeint/stepper/runge_kutta_dopri5.hpp>
#include <boost/numeric/odeint/stepper/bulirsch_stoer.hpp>

Eigen::VectorXd operator+ (const double s, const Eigen::VectorXd& v) {
    return (v.array() + s).matrix();
}

Eigen::VectorXd operator/ (const Eigen::VectorXd& v, const Eigen::VectorXd& w) {
    return (v.array() / w.array()).matrix();
}

namespace boost { namespace numeric { namespace odeint {
// create struct specialization for Eigen::VectorXd of
// the l-infinity norm used by controlled_runge_kutter
template<>
struct vector_space_norm_inf<Eigen::VectorXd>
{
    typedef double result_type; // = Eigen::VectorXd::Scalar
    double operator() (Eigen::VectorXd x) const
    {
        return x.lpNorm<Eigen::Infinity>();
    }
};
}}} // namespace boost::numeric::odeint

// Wrappers implementing IntegratorCore for boost::numeric::odeint Steppers (Controlled, Explicit and adams_bashforth_moulton)

namespace mio
{

template<
    size_t Steps,
    template<
        size_t InternalSteps = Steps,
        class State = Eigen::VectorXd,
        class Value = double,
        class Deriv = State,
        class Time = Value,
        class Algebra = boost::numeric::odeint::vector_space_algebra,
        class Operations = typename boost::numeric::odeint::operations_dispatcher<State>::operations_type,
        class Resizer = boost::numeric::odeint::never_resizer,
        class InitializingStepper = boost::numeric::odeint::runge_kutta4<State, Value, Deriv, Time, Algebra, Operations, Resizer>
    > class ABMStepper
>
class ABMStepperWrapper : public mio::IntegratorCore {
public:
    bool step(const mio::DerivFunction& f, Eigen::Ref<Eigen::VectorXd const> yt, double& t, double& dt,
              Eigen::Ref<Eigen::VectorXd> ytp1) const override {
        // reorder arguments of the DerivFunction f for the stepper
        std::function<void(const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, double t)> sys =
            [&](const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, double t){
                dxdt.resizeLike(x); // do_step calls sys with a vector of size 0 for some reason
                f(x, t, dxdt);
            };
        // copy y(t) to dxdt, since we use the scheme do_step(sys, inout, t, dt) with sys=f, inout=y(t) for
        // in-place computation - also, this form is shared by several (all?) steppers in boost
        Eigen::VectorXd dxdt = yt.eval();
        /* ABM has Steps internal states. according to the odeint guide on boost.org the stepper
         * should be called multiple times, but this loop gives incorrect results for the ode
        m_stepper.reset();
        for (size_t i = 0; i <= Steps; i++) {
            m_stepper.do_step(sys, dxdt, t, dt / Steps);
            //m_stepper.do_step(sys, dxdt, t, dt);
        }*/
        m_stepper.do_step(
            // reorder arguments of the DerivFunction f for the stepper
            [&](const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, double t){
                dxdt.resizeLike(x); // do_step calls sys with a vector of size 0 for some reason
                f(x, t, dxdt);
            },
            dxdt,
            t,
            dt
        );
        ytp1 = dxdt;
        t += dt;
        return true; // no step size adaption
    }
private:
    mutable ABMStepper<> m_stepper;
};

template<
    template<
        class State = Eigen::VectorXd,
        class Value = double,
        class Deriv = State,
        class Time = double,
        class Algebra = boost::numeric::odeint::vector_space_algebra,
        class Operations = typename boost::numeric::odeint::operations_dispatcher<State>::operations_type,
        class Resizer = boost::numeric::odeint::never_resizer
    > class ExplicitStepper
>
class ExplicitStepperWrapper : public mio::IntegratorCore {
public:
    bool step(const mio::DerivFunction& f, Eigen::Ref<Eigen::VectorXd const> yt, double& t, double& dt,
              Eigen::Ref<Eigen::VectorXd> ytp1) const override {
        // copy y(t) to dxdt, since we use the scheme do_step(sys, inout, t, dt) with sys=f, inout=y(t) for
        // in-place computation - also, this form is shared by several (all?) steppers in boost
        Eigen::VectorXd dxdt = yt.eval();
        m_stepper.do_step(
            // reorder arguments of the DerivFunction f for the stepper
            [&](const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, double t){
                dxdt.resizeLike(x); // do_step calls sys with a vector of size 0 for some reason
                f(x, t, dxdt);
            },
            dxdt,
            t,
            dt
        );
        ytp1 = dxdt;
        t += dt;
        return true; // no step size adaption
    }
private:
    mutable ExplicitStepper<> m_stepper;
};

template<
    template<
        class State = Eigen::VectorXd,
        class Value = double,
        class Deriv = State,
        class Time = double,
        class Algebra = boost::numeric::odeint::vector_space_algebra,
        class Operations = typename boost::numeric::odeint::operations_dispatcher<State>::operations_type,
        class Resizer = boost::numeric::odeint::never_resizer
    > class ControlledStepper
>
class ControlledStepperWrapper : public mio::IntegratorCore {
public:
    ControlledStepperWrapper(double abs_tol = 1e-6,
                             double rel_tol = 1e-6,
                             double dt_min=std::numeric_limits<double>::min(),
                             double dt_max=std::numeric_limits<double>::max()) :
        m_dt_min(dt_min),
        m_stepper(
            boost::numeric::odeint::default_error_checker<
                typename ControlledStepper<>::value_type,
                typename ControlledStepper<>::algebra_type,
                typename ControlledStepper<>::operations_type
            > (abs_tol, rel_tol),
            boost::numeric::odeint::default_step_adjuster<
                typename ControlledStepper<>::value_type,
                typename ControlledStepper<>::time_type
            > (dt_max)
        )
        // for more options see: boost/boost/numeric/odeint/stepper/controlled_runge_kutta.hpp
    {}
    bool step(const mio::DerivFunction& f, Eigen::Ref<Eigen::VectorXd const> yt, double& t, double& dt,
              Eigen::Ref<Eigen::VectorXd> ytp1) const override {
        // copy y(t) to dxdt, since we use the scheme try_step(sys, inout, t, dt) with sys=f, inout=y(t) for
        // in-place computation. This is similiar to do_step, but it updates t and dt
        Eigen::VectorXd dxdt = yt.eval();
        const double t_old = t; // t is updated by try_step on a successfull step
        while (t == t_old && dt > m_dt_min) {
            m_stepper.try_step(
                // reorder arguments of the DerivFunction f for the stepper
                [&](const Eigen::VectorXd& x, Eigen::VectorXd& dxdt, double t){
                    dxdt.resizeLike(x); // do_step calls sys with a vector of size 0 for some reason
                    f(x, t, dxdt);
                },
                dxdt,
                t,
                dt
            );
        }
        ytp1 = dxdt;
        return dt > m_dt_min;
    }
private:
    const double m_dt_min;
    mutable boost::numeric::odeint::controlled_runge_kutta<ControlledStepper<>> m_stepper;
};

} // namespace mio

#endif