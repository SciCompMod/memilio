/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#include "abm_helpers.h"
#include "actions.h"
#include "memilio/math/euler.h"
#include "memilio/math/adapt_rk.h"
#include "memilio/math/euler_maruyama.h"
#include "memilio/math/integrator.h"
#include "memilio/math/stepper_wrapper.h"
#include "memilio/utils/logging.h"
#include "utils.h"

#include "boost/numeric/odeint/stepper/modified_midpoint.hpp"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <vector>
#include <cmath>

void sin_deriv(Eigen::Ref<Eigen::VectorXd const> /*y*/, const double t, Eigen::Ref<Eigen::VectorXd> dydt)
{
    dydt[0] = std::cos(t);
}

template <class TestType>
class TestVerifyNumericalIntegrator : public testing::Test
{
protected:
    void SetUp() override
    {
        t0   = 0.;
        tmax = 2 * std::acos(-1); // 2PI
        err  = 0;
        mio::set_log_level(mio::LogLevel::off);
    }

    void TearDown() override
    {
        mio::set_log_level(mio::LogLevel::warn);
    }

public:
    std::vector<Eigen::VectorXd> y;
    std::vector<Eigen::VectorXd> sol;

    double t0;
    double tmax;
    size_t n;
    double dt;
    double err;
};

using ExplicitSteppersTestTypes =
    ::testing::Types<mio::EulerIntegratorCore<double>,
                     mio::ExplicitStepperWrapper<double, boost::numeric::odeint::modified_midpoint>,
                     mio::ExplicitStepperWrapper<double, boost::numeric::odeint::runge_kutta_cash_karp54>,
                     mio::ExplicitStepperWrapper<double, boost::numeric::odeint::runge_kutta_fehlberg78>>;

template <class T>
using TestVerifyExplicitNumericalIntegrator = TestVerifyNumericalIntegrator<::testing::Types<T>>;

TYPED_TEST_SUITE(TestVerifyExplicitNumericalIntegrator, ExplicitSteppersTestTypes);

TYPED_TEST(TestVerifyExplicitNumericalIntegrator, sine)
{
    this->n   = 1000;
    this->dt  = (this->tmax - this->t0) / this->n;
    this->y   = std::vector<Eigen::VectorXd>(this->n, Eigen::VectorXd::Constant(1, 0));
    this->sol = std::vector<Eigen::VectorXd>(this->n, Eigen::VectorXd::Constant(1, 0));

    this->sol[0][0]           = std::sin(0);
    this->sol[this->n - 1][0] = std::sin((this->n - 1) * this->dt);

    auto f = [](auto&& /*y*/, auto&& t, auto&& dydt) {
        dydt[0] = std::cos(t);
    };
    TypeParam stepper;

    auto t = this->t0;
    for (size_t i = 0; i < this->n - 1; i++) {
        this->sol[i + 1][0] = std::sin((i + 1) * this->dt);

        stepper.step(f, this->y[i], t, this->dt, this->y[i + 1]);

        this->err += std::pow(std::abs(this->y[i + 1][0] - this->sol[i + 1][0]), 2.0);
    }

    this->err = std::sqrt(this->err) / this->n;

    EXPECT_NEAR(this->err, 0.0, 1e-3);
}

using TestTypes = ::testing::Types<
    mio::RKIntegratorCore<double>,
    mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_cash_karp54>,
    // mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_dopri5>, // TODO: reenable once boost bug is fixed
    mio::ControlledStepperWrapper<double, boost::numeric::odeint::runge_kutta_fehlberg78>>;

TYPED_TEST_SUITE(TestVerifyNumericalIntegrator, TestTypes);

TYPED_TEST(TestVerifyNumericalIntegrator, sine)
{
    this->n   = 10;
    this->dt  = (this->tmax - this->t0) / this->n;
    this->y   = std::vector<Eigen::VectorXd>(this->n, Eigen::VectorXd::Constant(1, 0));
    this->sol = std::vector<Eigen::VectorXd>(this->n, Eigen::VectorXd::Constant(1, 0));

    TypeParam rk;
    rk.set_abs_tolerance(1e-7);
    rk.set_rel_tolerance(1e-7);
    rk.set_dt_min(1e-3);
    rk.set_dt_max(1.0);

    this->sol[0][0] = std::sin(0);

    size_t i      = 0;
    double t_eval = this->t0;
    // printf("\n t: %.8f\t sol %.8f\t rkf %.8f", t, sol[0][0], y[0][0]);
    while (t_eval - this->tmax < 1e-10) {

        if (i + 1 >= this->sol.size()) {
            this->sol.push_back(Eigen::VectorXd::Constant(1, 0));
            this->y.push_back(Eigen::VectorXd::Constant(1, 0));
        }

        rk.step(&sin_deriv, this->y[i], t_eval, this->dt, this->y[i + 1]); //

        this->sol[i + 1][0] = std::sin(t_eval);

        // printf("\n t: %.8f (dt %.8f)\t sol %.8f\t rkf %.8f", t_eval, dt, sol[i + 1][0], y[i + 1][0]);
        // printf("\n approx: %.4e, sol: %.4e, error %.4e", y[i+1][0], sol[i+1][0], err);

        this->err += std::pow(std::abs(this->y[i + 1][0] - this->sol[i + 1][0]), 2.0);
        i++;
    }

    this->n = i;

    this->err = std::sqrt(this->err) / this->n;

    EXPECT_NEAR(this->err, 0.0, 1e-7);
}

TYPED_TEST(TestVerifyNumericalIntegrator, adaptiveStepSizing)
{
    // this test checks all requirements on adaptive step sizing for IntegratorCore::step
    // set fixed tolarances to control when a integration step fails
    const double tol    = 1;
    const double dt_min = 1, dt_max = 2 * dt_min;

    this->y   = {Eigen::VectorXd::Zero(1)};
    this->sol = {Eigen::VectorXd::Zero(1)};

    TypeParam integrator;
    integrator.set_abs_tolerance(tol);
    integrator.set_rel_tolerance(tol);
    integrator.set_dt_min(dt_min);
    integrator.set_dt_max(dt_max);

    double t_eval;
    bool step_okay;

    // this deriv function is supposed to (not guaranteed to!) break any integrator
    double c        = 1;
    auto deriv_fail = [&c, tol](const auto&&, auto&&, auto&& dxds) {
        c /= -10 * tol;
        dxds.array() = 1 / c; // increasing oscillation with each evaluation (indep. of t and dt)
    };
    // this deriv function is easily integrable
    auto deriv_success = [](const auto&&, auto&&, auto&& dxds) {
        dxds.array() = 1; // const f
    };

    // check that on a failed step, dt does decrease, but not below dt_min

    t_eval    = 0;
    this->dt  = dt_max;
    step_okay = integrator.step(deriv_fail, this->y[0], t_eval, this->dt, this->sol[0]);
    c         = 1; // reset deriv_fail

    EXPECT_EQ(step_okay, false); // step sizing should fail
    EXPECT_EQ(this->dt, dt_min); // new step size should fall back to dt_min
    EXPECT_EQ(t_eval, dt_min); // a step must be made with no less than dt_min

    // check that on a successful step, dt increases from dt_min

    t_eval    = 0;
    this->dt  = dt_min;
    step_okay = integrator.step(deriv_success, this->y[0], t_eval, this->dt, this->sol[0]);

    EXPECT_EQ(step_okay, true); // step sizing must be okay
    EXPECT_GE(this->dt, dt_min); // new step size may be larger
    EXPECT_LE(this->dt, dt_max); // but not too large
    EXPECT_EQ(t_eval, dt_min); // used step size should still be dt_min
    EXPECT_EQ(this->sol[0][0], t_eval); // check that the integration step matches the time step

    // check that on a successful step, dt does not increase from dt_max

    t_eval    = 0;
    this->dt  = dt_max;
    step_okay = integrator.step(deriv_success, this->y[0], t_eval, this->dt, this->sol[0]);

    EXPECT_EQ(step_okay, true); // step sizing must be okay
    EXPECT_EQ(this->dt, dt_max); // new step size must not increase from dt_max
    EXPECT_EQ(t_eval, dt_max); // used step size should be dt_max
    EXPECT_EQ(this->sol[0][0], t_eval); // check that the integration step matches the time step

    // check that the integrator does not make time steps for dt not in [dt_min, dt_max]

    t_eval    = 0;
    this->dt  = 0.5 * dt_min;
    step_okay = integrator.step(deriv_success, this->y[0], t_eval, this->dt, this->sol[0]);

    EXPECT_EQ(step_okay, true); // step sizing must be okay
    EXPECT_GE(this->dt, dt_min); // new step size must be bounded from below
    EXPECT_LE(this->dt, dt_max); // and above
    EXPECT_EQ(t_eval, dt_min); // step size should fall back to dt_min
    EXPECT_EQ(this->sol[0][0], t_eval); // check that the integration step matches the time step

    t_eval    = 0;
    this->dt  = 2.0 * dt_max;
    step_okay = integrator.step(deriv_success, this->y[0], t_eval, this->dt, this->sol[0]);

    EXPECT_EQ(step_okay, true); // step sizing must be okay
    EXPECT_GE(this->dt, dt_min); // new step size must be bounded from below
    EXPECT_LE(this->dt, dt_max); // and above
    EXPECT_EQ(t_eval, dt_max); // step size should fall back to dt_max
    EXPECT_EQ(this->sol[0][0], t_eval); // check that the integration step matches the time step
}

auto DoStep()
{
    return testing::DoAll(testing::WithArgs<2, 3>(AddAssign()), testing::WithArgs<4, 1>(AssignUnsafe()),
                          testing::Return(true));
}

class MockIntegratorCore : public mio::OdeIntegratorCore<double>
{
public:
    MockIntegratorCore()
        : MockIntegratorCore(0.0, 0.0)
    {
    }

    MockIntegratorCore(double dt_min, double dt_max)
        : mio::OdeIntegratorCore<double>(dt_min, dt_max)
    {
        ON_CALL(*this, step).WillByDefault(DoStep());
    }
    MOCK_METHOD(bool, step,
                (const mio::DerivFunction<double>& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
                 Eigen::Ref<Eigen::VectorXd> ytp1),
                (const));
    
    std::unique_ptr<mio::OdeIntegratorCore<double>> clone() const override 
    {
        throw std::runtime_error("MockIntegratorCore clone() called unexpectedly");
    }
};

TEST(TestOdeIntegrator, integratorDoesTheRightNumberOfSteps)
{
    using testing::_;
    auto mock_core = std::make_unique<testing::StrictMock<MockIntegratorCore>>();
    EXPECT_CALL(*mock_core, step).Times(100);

    auto f          = [](auto&&, auto&&, auto&&) {};
    auto integrator = mio::OdeIntegrator<double>(std::move(mock_core));
    mio::TimeSeries<double> result(0, Eigen::VectorXd::Constant(1, 0.0));
    double dt = 1e-2;
    integrator.advance(f, 1, dt, result);
    EXPECT_EQ(result.get_num_time_points(), 101);
}

TEST(TestOdeIntegrator, integratorStopsAtTMax)
{
    auto f = [](auto&&, auto&&, auto&&) {};
    mio::TimeSeries<double> result(0, Eigen::VectorXd::Constant(1, 0.0));
    double dt       = 0.137;
    auto integrator = mio::OdeIntegrator<double>(std::make_unique<testing::NiceMock<MockIntegratorCore>>());
    integrator.advance(f, 2.34, dt, result);
    EXPECT_DOUBLE_EQ(result.get_last_time(), 2.34);
}

auto DoStepAndIncreaseStepsize(double new_dt)
{
    return testing::DoAll(testing::WithArgs<2, 3>(AddAssign()), testing::WithArgs<4, 1>(AssignUnsafe()),
                          testing::SetArgReferee<3>(new_dt), testing::Return(true));
}

auto DoStepAndReduceStepsize(double new_dt)
{
    return testing::DoAll(testing::WithArgs<2>(AddAssign(new_dt)), testing::WithArgs<4, 1>(AssignUnsafe()),
                          testing::SetArgReferee<3>(new_dt), testing::Return(true));
}

TEST(TestOdeIntegrator, integratorUpdatesStepsize)
{
    using testing::_;
    using testing::Eq;

    auto mock_core = std::make_unique<testing::StrictMock<MockIntegratorCore>>();

    {
        testing::InSequence seq;

        EXPECT_CALL(*mock_core, step(_, _, _, Eq(1), _)).Times(1).WillOnce(DoStepAndIncreaseStepsize(2.));
        EXPECT_CALL(*mock_core, step(_, _, _, Eq(2), _)).Times(1).WillOnce(DoStepAndIncreaseStepsize(4.));
        EXPECT_CALL(*mock_core, step(_, _, _, Eq(4), _)).Times(1).WillOnce(DoStepAndIncreaseStepsize(8.));
        //last step of first advance call
        EXPECT_CALL(*mock_core, step(_, _, _, Eq(3), _))
            .Times(1)
            .WillRepeatedly(DoStepAndIncreaseStepsize(6.)); //this stepsize should not be stored!
        //continue on the second advance call with correct stepsize
        EXPECT_CALL(*mock_core, step(_, _, _, Eq(8), _)).Times(1).WillOnce(DoStepAndReduceStepsize(3.5));
        EXPECT_CALL(*mock_core, step(_, _, _, Eq(3.5), _)).Times(1).WillOnce(DoStepAndIncreaseStepsize(8));
        //last step of second advance call
        EXPECT_CALL(*mock_core, step(_, _, _, Eq(6), _)).Times(1);
    }

    auto f = [](auto&&, auto&&, auto&&) {};
    mio::TimeSeries<double> result(0, Eigen::VectorXd::Constant(1, 0.0));
    double dt       = 1.0;
    auto integrator = mio::OdeIntegrator<double>(std::move(mock_core));
    integrator.advance(f, 10.0, dt, result);
    integrator.advance(f, 23.0, dt, result);
}

auto DoStepAndIncreaseY(const Eigen::VectorXd& dy)
{
    return testing::DoAll(testing::WithArgs<2, 3>(AddAssign()), testing::WithArgs<4, 1>(AssignUnsafe()),
                          testing::WithArgs<4>(AddAssignUnsafe(dy)), testing::Return(true));
}

TEST(TestOdeIntegrator, integratorContinuesAtLastState)
{
    using testing::_;
    using testing::Eq;

    double dt0      = 0.25;
    double dt       = dt0;
    auto dy         = Eigen::VectorXd::Constant(1, 1);
    auto y0         = Eigen::VectorXd::Constant(1, 0);
    auto mock_core  = std::make_unique<testing::StrictMock<MockIntegratorCore>>();

    {
        testing::InSequence seq;
        EXPECT_CALL(*mock_core, step(_, Eq(y0), _, _, _)).WillOnce(DoStepAndIncreaseY(dy));
        EXPECT_CALL(*mock_core, step(_, Eq(y0 + dy), _, _, _)).WillOnce(DoStepAndIncreaseY(dy));
        EXPECT_CALL(*mock_core, step(_, Eq(y0 + 2 * dy), _, _, _)).WillOnce(DoStepAndIncreaseY(dy));
        EXPECT_CALL(*mock_core, step(_, Eq(y0 + 3 * dy), _, _, _)).WillOnce(DoStepAndIncreaseY(dy));
        EXPECT_CALL(*mock_core, step(_, Eq(y0 + 4 * dy), _, _, _)).WillOnce(DoStepAndIncreaseY(dy));
    }

    auto f          = [](auto&&, auto&&, auto&&) {};
    auto integrator = mio::OdeIntegrator<double>(std::move(mock_core));
    mio::TimeSeries<double> result(0, y0);
    integrator.advance(f, 4 * dt0, dt, result);
    integrator.advance(f, 5 * dt0, dt, result);
}

TEST(TestOdeIntegrator, integratorForcesLastStepSize)
{
    using testing::_;
    using testing::Eq;

    const double dt_min  = 0.7;
    double dt            = 0.5; // this is on purpose < dt_min
    const double t_max   = 3.0; // must not be an integer multiple of dt_min
    auto mock_core       = std::make_unique<testing::StrictMock<MockIntegratorCore>>(dt_min, t_max);
    auto f               = [](auto&&, auto&&, auto&&) {};
    auto step_fct        = [&mock = *mock_core](auto&&, auto&&, auto& t_, auto& dt_, auto&&) {
        dt_ = std::max(dt_, mock.get_dt_min());
        t_ += dt_;
        return true;
    };


    const auto num_calls = Eigen::Index(t_max / dt_min) + 1;
    EXPECT_CALL(*mock_core, step).Times(num_calls).WillRepeatedly(testing::Invoke(step_fct));

    // run a mock integration to examine whether only the last step is forced
    mio::TimeSeries<double> mock_result(0, Eigen::VectorXd::Constant(1, 0));
    auto integrator      = mio::OdeIntegrator<double>(std::move(mock_core));
    integrator.advance(f, t_max, dt, mock_result);
    EXPECT_EQ(mock_result.get_num_time_points(), num_calls + 1);
    for (Eigen::Index i = 0; i < mock_result.get_num_time_points(); i++) {
        EXPECT_DOUBLE_EQ(mock_result.get_time(i), std::min(i * dt_min, t_max));
    }
}

TEST(TestStochasticIntegrator, EulerMaruyamaIntegratorCore)
{
    using X = Eigen::Vector2d;

    mio::EulerMaruyamaIntegratorCore<double> integrator;
    const double ref_val   = std::sin(1);
    const double ref_noise = std::sqrt(ref_val);
    double white_noise     = 1;
    double t, dt = 1;
    mio::RedirectLogger logger;

    const auto f = [](Eigen::Ref<const X> _, double s, Eigen::Ref<X> y) {
        mio::unused(_);
        y[0] = std::sin(s);
        y[1] = 0;
    };
    const auto noise = [&f, &white_noise](Eigen::Ref<const X> x, double s, Eigen::Ref<X> y) {
        f(x, s, y);
        y = white_noise * y.array().sqrt();
    };

    logger.capture();

    X x = {0, 1}, y; // feasable setup, with x[1] as rescale target
    t = 1, white_noise = 1;
    integrator.step(f, noise, x, t, dt, y); // no negative value, no rescaling
    EXPECT_NEAR(y[0], ref_val + ref_noise, mio::Limits<double>::zero_tolerance());
    EXPECT_NEAR(y[1], 1, mio::Limits<double>::zero_tolerance());

    t = 1, white_noise = -1;
    integrator.step(f, noise, x, t, dt, y); // rescale succeeds
    EXPECT_NEAR(y[0], 0, mio::Limits<double>::zero_tolerance());
    EXPECT_NEAR(y[1], 1 + (ref_val - ref_noise), mio::Limits<double>::zero_tolerance());

    EXPECT_TRUE(logger.read().empty()); // no log messages yet

    x[1] = 0; // take away rescale target
    t = 1, white_noise = -1;
    integrator.step(f, noise, x, t, dt, y); // no positive values, rescale fails
    EXPECT_GT(y[0], 0);
    EXPECT_NEAR(y[0], y[1], mio::Limits<double>::zero_tolerance());
    EXPECT_THAT(logger.read(), testing::HasSubstr("[error] Failed to rescale values"));

    t = 0, white_noise = -1;
    integrator.step(f, noise, x, t, dt, y); // everything is 0, rescale does nothing
    EXPECT_EQ(y[0], 0);
    EXPECT_EQ(y[1], 0);
    EXPECT_TRUE(logger.read().empty()); // no log messages again

    logger.release();
}
