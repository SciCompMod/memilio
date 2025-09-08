/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Ralf Hannemann-Tamas, Lena Ploetzke
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
#include "ode_seair/model.h"
#include "ode_seair/infection_state.h"
#include "ode_seair/parameters.h"
#include "memilio/compartments/simulation.h"

#include "memilio/utils/logging.h"

#include "ad/ad_wrapper.hpp"
#include "boost/numeric/odeint.hpp"
#include <gtest/gtest.h>

template <typename FP>
FP square(FP x)
{
    return x * x;
}

// Test that the "ad" library can be used as expected in the software memilio.
TEST(Testad, ad_square)
{
    std::vector<double> x_values = {3., 5., 10.5};
    for (const double& x : x_values) {
        // Forward AD.
        ad::gt1s<double>::type t1_x; // tangent-linear AD type
        ad::value(t1_x)             = x; // set value (this corresponds to the non-AD part)
        ad::derivative(t1_x)        = 1.0; // set tangent-linear derivative seed
        ad::gt1s<double>::type t1_y = square(t1_x);
        EXPECT_NEAR(ad::value(t1_x) * ad::value(t1_x), ad::value(t1_y), 1e-8);
        EXPECT_NEAR(2. * ad::value(t1_x), ad::derivative(t1_y), 1e-8);

        // Reverse AD.
        // Create tape (allocation).
        if (!ad::ga1s<double>::global_tape)
            ad::ga1s<double>::global_tape = ad::ga1s<double>::tape_t::create();
        // Clear tape.
        ad::ga1s<double>::global_tape->reset();

        ad::ga1s<double>::type a1_x; // Reverse-mode AD type.
        ad::value(a1_x)      = x; // Set value (this corresponds to the non-AD part)
        ad::derivative(a1_x) = 0.0; // For reverse-mode the seed has to be set to zero!
        ad::ga1s<double>::global_tape->register_variable(a1_x); // Register input variable.
        ad::ga1s<double>::type a1_y = square(a1_x);

        EXPECT_NEAR(ad::value(a1_x) * ad::value(a1_x), ad::value(a1_y), 1e-8);

        ad::ga1s<double>::global_tape->register_output_variable(a1_y);
        ad::derivative(a1_y) = 1.0;
        ad::ga1s<double>::global_tape
            ->interpret_adjoint(); // Compute reverse-mode derivatives by evaluating the tape backwards.
        // Access reverse-derivatives in x variable.
        EXPECT_NEAR(2 * ad::value(a1_x), ad::derivative(a1_x), 1e-8);
        ad::ga1s<double>::tape_t::remove(ad::ga1s<double>::global_tape); // Deallocate tape.
    }
}

// Test that the "ad" library can be used as expected with a more complex example.
// This test ensures that boost::numeric::odeint::runge_kutta_cash_karp54 can be fully
// algorithmically diffentiated using the algorithmic differentiation (AD) data types of ad/ad.hpp.

// Define the rhs of the ODE x' = f(x), that should be solved using odeint.
template <typename value_type, typename time_type>
void harmonic_oscillator(const std::array<value_type, 2>& x, std::array<value_type, 2>& dxdt, const time_type /* t */)
{
    const double damping = 0.15;
    dxdt[0]              = x[1];
    dxdt[1]              = -x[0] - damping * x[1];
}

TEST(Testad, ad_odeint)
{
    using ad_forward_t = typename ad::gt1s<double>::type; // AD data type for scalar forward mode.
    // The type of container used to hold the state vector.
    using value_type = ad_forward_t;
    using time_type  = value_type;
    using state_type = std::array<value_type, 2>; // 2-dimensional vector

    using error_stepper_type =
        boost::numeric::odeint::runge_kutta_cash_karp54<state_type, value_type, state_type, time_type>;

    state_type x;
    x[0]                 = 1.0; // Start at x=1.0, y=0.0.
    x[1]                 = 0.0;
    ad::derivative(x[0]) = 1.0; // Compute derivative with respect to x[0] (scalar tangent-linear mode).

    auto t0    = time_type(0.0);
    auto t_end = time_type(10.0);
    auto dt    = time_type(0.01);

    const double abs_tol = 1e-6; // Absolute tolerance for error-controlled integration.
    const double rel_tol = 1e-6; // Relative tolerance for error-controlled integration.

    boost::numeric::odeint::integrate_adaptive(
        boost::numeric::odeint::make_controlled<error_stepper_type>(abs_tol, rel_tol),
        harmonic_oscillator<value_type, time_type>, x, t0, t_end, dt);

    // We want to compare AD derivatives with difference quotient.
    // To this end, we simulate again with a small perturbation of the initial value of x[0].
    const double h = 1e-3;
    state_type x_compare;
    x_compare[0] = 1.0 + h;
    x_compare[1] = 0.0;
    boost::numeric::odeint::integrate_adaptive(
        boost::numeric::odeint::make_controlled<error_stepper_type>(abs_tol, rel_tol),
        harmonic_oscillator<value_type, time_type>, x_compare, t0, t_end, dt);

    // Compare AD derivatives with a difference quotient.
    EXPECT_NEAR(ad::derivative(x[0]), (ad::value(x_compare[0]) - ad::value(x[0])) / h, 1e-4);
    EXPECT_NEAR(ad::derivative(x[1]), (ad::value(x_compare[1]) - ad::value(x[1])) / h, 1e-4);
}

// Check that ad correctly takes the derivative of an ODE model (according to a parameter).
TEST(Testad, ad_seair)
{
    using FP = typename ad::gt1s<double>::type; // AD data type for scalar forward mode.

    FP t0   = 0.;
    FP tmax = 0.1;
    FP dt   = 0.1;

    // Define model with AD data type.
    mio::oseair::Model<FP> admodel;

    // Set initial values.
    admodel.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}]  = 450.;
    admodel.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}]      = 100.;
    admodel.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}] = 200.;
    admodel.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}]     = 50.;
    admodel.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}]    = 100.;
    admodel.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Dead)}]         = 100.;

    // Compute derivative with respect to the TestingRate (scalar tangent-linear mode).
    ad::derivative(admodel.parameters.get<mio::oseair::TestingRate<FP>>()) = 1.0;
    ad::value(admodel.parameters.get<mio::oseair::TestingRate<FP>>())      = 0.2;

    auto adresult = mio::simulate<FP, mio::oseair::Model<FP>>(t0, tmax, dt, admodel);

    // We want to compare AD derivatives with a difference quotient.
    // Therefore, define a model with double data type with a small pertubation of the parameter TestingRate.
    const double h = 1e-4;
    mio::oseair::Model<double> model;

    // Set same initial values as above.
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Susceptible)}]  = 450.;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Exposed)}]      = 100.;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Asymptomatic)}] = 200.;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Infected)}]     = 50.;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Recovered)}]    = 100.;
    model.populations[{mio::Index<mio::oseair::InfectionState>(mio::oseair::InfectionState::Dead)}]         = 100.;

    // Small pertubation of the parameter TestingRate.
    model.parameters.get<mio::oseair::TestingRate<double>>() =
        ad::value(admodel.parameters.get<mio::oseair::TestingRate<FP>>()) + h;

    auto result =
        mio::simulate<double, mio::oseair::Model<double>>(ad::value(t0), ad::value(tmax), ad::value(dt), model);

    // Check that the derivative calculated with ad is close to the result obtained with a difference quotient.
    // As we have an adaptive method, we can also test the influence on the simulation time.
    EXPECT_NEAR(ad::derivative(adresult.get_last_time()),
                (result.get_last_time() - ad::value(adresult.get_last_time())) / h, 1e-3);
    // Derivative of the compartment sizes.
    for (int i = 0; i < (int)mio::oseair::InfectionState::Count; i++) {
        EXPECT_NEAR(ad::derivative(adresult.get_last_value()[i]),
                    (result.get_last_value()[i] - ad::value(adresult.get_last_value()[i])) / h, 1e-3);
    }
}
