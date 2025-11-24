/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Ralf Hannemann-Tamas
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

#include "memilio/ad/ad.h"
#include "boost/numeric/odeint.hpp"
#include <array>

// This program shows that  boost::numeric::odeint::runge_kutta_cash_karp54 can be fully
// algorithmically diffentiated using the algorithmic differentiation (AD) data types of ad/ad.hpp.

using ad_forward_type = typename ad::gt1s<double>::type; // AD data type for scalar forward mode

using value_type = ad_forward_type;
using time_type  = value_type;

// The type of container used to hold the state vector
const size_t dimension = 2; // dimension of ODE
using state_type       = std::array<value_type, dimension>; // 2-dimensional vector

double damping = 0.15;

/** The rhs of x' = f(x). This ODE describes a damped harmonic oscillator. */
void harmonic_oscillator(const state_type& x, state_type& dxdt, const time_type /* t */)
{
    dxdt[0] = x[1];
    dxdt[1] = -x[0] - damping * x[1];
}

// We use a an adaptive step size control
using error_stepper_type =
    boost::numeric::odeint::runge_kutta_cash_karp54<state_type, value_type, state_type, value_type>;
using controlled_stepper_type = boost::numeric::odeint::controlled_runge_kutta<error_stepper_type>;

using boost::numeric::odeint::make_controlled;

int main()
{
    state_type x;
    x[0]                 = 1.0; // start at x=1.0, p=0.0
    x[1]                 = 0.0;
    ad::derivative(x[0]) = 1.0; // compute derivative with respect to x[0] (scalar tangent-linear mode)
    ad::derivative(x[1]) = 0.0;

    auto t0    = time_type(0.0); // initial time
    auto t_end = time_type(10.0); // stop time
    auto dt    = time_type(0.01); // hint for initial step size

    const double abs_tol = 1e-6; // absolute tolerance for error-controlled integration
    const double rel_tol = 1e-6; // relative tolerance for error-controlled integration

    controlled_stepper_type controlled_stepper;
    boost::numeric::odeint::integrate_adaptive(make_controlled<error_stepper_type>(abs_tol, rel_tol),
                                               harmonic_oscillator, x, t0, t_end, dt);
    std::cout << "ad derivative of x is (" << ad::derivative(x[0]) << ", " << ad::derivative(x[1]) << ")\n";

    // We want to compare AD derivatives with difference quotient
    // To this end, we simulate again with a perturbation of the initial value of x[0]
    const double h                  = 1e-3; // pertubation for finite differences
    std::array<double, dimension> y = {ad::value(x[0]), ad::value(x[1])};
    x[0]                            = 1.0 + h; // add perturbation to initial value of x[0]
    x[1]                            = 0.0;

    // integrate perturbed system
    boost::numeric::odeint::integrate_adaptive(make_controlled<error_stepper_type>(abs_tol, rel_tol),
                                               harmonic_oscillator, x, t0, t_end, dt);

    std::cout << "finite differences derivative of x is (" << (ad::value(x[0]) - y[0]) / h << ", "
              << (ad::value(x[1]) - y[1]) / h << ")\n";

    return 0;
}
