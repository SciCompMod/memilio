/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "memilio/math/euler.h"
#include "memilio/math/adapt_rk.h"
#include "memilio/math/stepper_wrapper.h"
#include <actions.h>

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <string>
#include <vector>
#include <fstream>
#include <ios>
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

using EulerTest = TestVerifyNumericalIntegrator<::testing::Types<mio::EulerIntegratorCore>>;
TEST_F(EulerTest, euler_sine)
{
    n   = 1000;
    dt  = (tmax - t0) / n;
    y   = std::vector<Eigen::VectorXd>(n, Eigen::VectorXd::Constant(1, 0));
    sol = std::vector<Eigen::VectorXd>(n, Eigen::VectorXd::Constant(1, 0));

    sol[0][0]     = std::sin(0);
    sol[n - 1][0] = std::sin((n - 1) * dt);

    auto f = [](auto&& /*y*/, auto&& t, auto&& dydt) {
        dydt[0] = std::cos(t);
    };
    mio::EulerIntegratorCore euler;

    auto t = t0;
    for (size_t i = 0; i < n - 1; i++) {
        sol[i + 1][0] = std::sin((i + 1) * dt);

        euler.step(f, y[i], t, dt, y[i + 1]);

        // printf("\n %.8f\t %.8f ", y[i + 1][0], sol[i + 1][0]);

        err += std::pow(std::abs(y[i + 1][0] - sol[i + 1][0]), 2.0);
    }

    err = std::sqrt(err) / n;

    EXPECT_NEAR(err, 0.0, 1e-3);
}

using TestTypes = ::testing::Types<mio::RKFastIntegratorCore, mio::RKIntegratorCore, mio::VRKFastIntegratorCore,
                                   mio::VRKFastIntegratorCore, mio::VRKOptIntegratorCore,
                                   mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>,
                                   mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_dopri5>,
                                   mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_fehlberg78>>;

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

auto DoStep()
{
    return testing::DoAll(testing::WithArgs<2, 3>(AddAssign()), testing::WithArgs<4, 1>(AssignUnsafe()),
                          testing::Return(true));
}

class MockIntegratorCore : public mio::IntegratorCore
{
public:
    MockIntegratorCore()
    {
        ON_CALL(*this, step).WillByDefault(DoStep());
    }
    MOCK_METHOD(bool, step,
                (const mio::DerivFunction& f, Eigen::Ref<const Eigen::VectorXd> yt, double& t, double& dt,
                 Eigen::Ref<Eigen::VectorXd> ytp1),
                (const));
};

TEST(TestOdeIntegrator, integratorDoesTheRightNumberOfSteps)
{
    using testing::_;
    auto mock_core = std::make_shared<testing::StrictMock<MockIntegratorCore>>();
    EXPECT_CALL(*mock_core, step).Times(100);

    auto f          = [](auto&&, auto&&, auto&&) {};
    auto integrator = mio::OdeIntegrator(f, 0, Eigen::VectorXd::Constant(1, 0.0), 1e-2, mock_core);
    integrator.advance(1);
    EXPECT_EQ(integrator.get_result().get_num_time_points(), 101);
}

TEST(TestOdeIntegrator, integratorStopsAtTMax)
{
    auto f          = [](auto&&, auto&&, auto&&) {};
    auto integrator = mio::OdeIntegrator(f, 0, Eigen::VectorXd::Constant(1, 0.0), 0.137,
                                         std::make_shared<testing::NiceMock<MockIntegratorCore>>());
    integrator.advance(2.34);
    EXPECT_DOUBLE_EQ(integrator.get_result().get_last_time(), 2.34);
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

    auto mock_core = std::make_shared<testing::StrictMock<MockIntegratorCore>>();

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

    auto f          = [](auto&&, auto&&, auto&&) {};
    auto integrator = mio::OdeIntegrator(f, 0, Eigen::VectorXd::Constant(1, 0), 1.0, mock_core);
    integrator.advance(10.0);
    integrator.advance(23.0);
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

    double dt       = 0.25;
    auto dy         = Eigen::VectorXd::Constant(1, 1);
    auto y0         = Eigen::VectorXd::Constant(1, 0);
    auto mock_core  = std::make_shared<testing::StrictMock<MockIntegratorCore>>();
    auto f          = [](auto&&, auto&&, auto&&) {};
    auto integrator = mio::OdeIntegrator(f, 0, y0, dt, mock_core);

    {
        testing::InSequence seq;
        EXPECT_CALL(*mock_core, step(_, Eq(y0), _, _, _)).WillOnce(DoStepAndIncreaseY(dy));
        EXPECT_CALL(*mock_core, step(_, Eq(y0 + dy), _, _, _)).WillOnce(DoStepAndIncreaseY(dy));
        EXPECT_CALL(*mock_core, step(_, Eq(y0 + 2 * dy), _, _, _)).WillOnce(DoStepAndIncreaseY(dy));
        EXPECT_CALL(*mock_core, step(_, Eq(y0 + 3 * dy), _, _, _)).WillOnce(DoStepAndIncreaseY(dy));
        EXPECT_CALL(*mock_core, step(_, Eq(y0 + 4 * dy), _, _, _)).WillOnce(DoStepAndIncreaseY(dy));
    }

    integrator.advance(4 * dt);
    integrator.advance(5 * dt);
}
