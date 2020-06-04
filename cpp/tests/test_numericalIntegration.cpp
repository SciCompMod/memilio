#include <gtest_helpers.h>

#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <epidemiology/euler.h>
#include <epidemiology/adapt_rk.h>

#include <string>
#include <vector>
#include <fstream>
#include <ios>
#include <cmath>

void sin_deriv(Eigen::VectorXd const& y, const double t, Eigen::VectorXd& dydt)
{
    dydt[0] = std::cos(t);
}

class TestVerifyNumericalIntegrator : public testing::Test
{
protected:
    void SetUp() override
    {
        t    = 0.;
        tmax = 2 * std::acos(-1); // 2PI
        err  = 0;
    }

public:
    std::vector<Eigen::VectorXd> y;
    std::vector<Eigen::VectorXd> sol;

    double t;
    double tmax;
    size_t n;
    double dt;
    double err;
};

TEST_F(TestVerifyNumericalIntegrator, euler_sine)
{
    n   = 1000;
    dt  = (tmax - t) / n;
    y   = std::vector<Eigen::VectorXd>(n, Eigen::VectorXd::Constant(1, 0));
    sol = std::vector<Eigen::VectorXd>(n, Eigen::VectorXd::Constant(1, 0));

    sol[0][0]     = std::sin(0);
    sol[n - 1][0] = std::sin((n - 1) * dt);

    auto f = [](Eigen::VectorXd const& y, const double t, Eigen::VectorXd& dydt) {
        dydt[0] = std::cos(t);
    };
    epi::EulerIntegratorCore euler;

    for (size_t i = 0; i < n - 1; i++) {
        sol[i + 1][0] = std::sin((i + 1) * dt);

        euler.step(f, y[i], t, dt, y[i + 1]);

        // printf("\n %.8f\t %.8f ", y[i + 1][0], sol[i + 1][0]);

        err += std::pow(std::abs(y[i + 1][0] - sol[i + 1][0]), 2.0);
    }

    err = std::sqrt(err) / n;

    EXPECT_NEAR(err, 0.0, 1e-3);
}

TEST_F(TestVerifyNumericalIntegrator, runge_kutta_fehlberg45_sine)
{

    n   = 10;
    dt  = (tmax - t) / n;
    y   = std::vector<Eigen::VectorXd>(n, Eigen::VectorXd::Constant(1, 0));
    sol = std::vector<Eigen::VectorXd>(n, Eigen::VectorXd::Constant(1, 0));

    epi::RKIntegratorCore rkf45(1e-3, 1.0);
    rkf45.set_abs_tolerance(1e-7);
    rkf45.set_rel_tolerance(1e-7);

    sol[0][0] = std::sin(0);

    size_t i      = 0;
    double t_eval = t;
    // printf("\n t: %.8f\t sol %.8f\t rkf %.8f", t, sol[0][0], y[0][0]);
    while (t_eval - tmax < 1e-10) {

        if (i + 1 >= sol.size()) {
            sol.push_back(Eigen::VectorXd::Constant(1, 0));
            y.push_back(Eigen::VectorXd::Constant(1, 0));
        }

        double dt_old = dt;

        rkf45.step(&sin_deriv, y[i], t_eval, dt, y[i + 1]); //

        sol[i + 1][0] = std::sin(t_eval);

        // printf("\n t: %.8f (dt %.8f)\t sol %.8f\t rkf %.8f", t_eval, dt, sol[i + 1][0], y[i + 1][0]);
        // printf("\n approx: %.4e, sol: %.4e, error %.4e", y[i+1][0], sol[i+1][0], err);

        err += std::pow(std::abs(y[i + 1][0] - sol[i + 1][0]), 2.0);
        i++;
    }

    n = i;

    err = std::sqrt(err) / n;

    EXPECT_NEAR(err, 0.0, 1e-7);
}

auto DoStep()
{
    return testing::DoAll(
        testing::WithArgs<2, 3>(AddAssign()),
        testing::WithArgs<4, 1>(Assign()),
        testing::Return(true));
}

class MockIntegratorCore : public epi::IntegratorCore
{
public:
    MockIntegratorCore()
    {
        ON_CALL(*this, step).WillByDefault(DoStep());
    }
    MOCK_METHOD(bool, step, (const epi::DerivFunction& f, const Eigen::VectorXd& yt, double& t, double& dt, Eigen::VectorXd& ytp1), (const));
};

TEST(TestOdeIntegrator, integratorDoesTheRightNumberOfSteps)
{
    using testing::_;
    auto mock_core = std::make_shared<MockIntegratorCore>();
    EXPECT_CALL(*mock_core, step).Times(100);

    auto f = [](const auto& y, auto t, auto& dydt) {};
    auto integrator = epi::OdeIntegrator(f, 0, Eigen::VectorXd::Constant(1, 0.0), 1e-2, mock_core);
    integrator.advance(1);
    EXPECT_EQ(integrator.get_t().size(), 101);
    EXPECT_EQ(integrator.get_y().size(), 101);
}

TEST(TestOdeIntegrator, integratorStopsAtTMax)
{
    auto f          = [](const auto& y, auto t, auto& dydt) {};
    auto integrator = epi::OdeIntegrator(f, 0, Eigen::VectorXd::Constant(1, 0.0), 0.137,
                                         std::make_shared<testing::NiceMock<MockIntegratorCore>>());
    integrator.advance(2.34);
    EXPECT_DOUBLE_EQ(integrator.get_t().back(), 2.34);
}

auto DoStepAndSetStepsize(double new_dt)
{
    return testing::DoAll(
        testing::WithArgs<2, 3>(AddAssign()),
        testing::WithArgs<4, 1>(Assign()),
        testing::SetArgReferee<3>(new_dt),
        testing::Return(true));
}

TEST(TestOdeIntegrator, integratorUpdatesStepsize)
{
    using testing::_;
    using testing::Eq;

    auto mock_core = std::make_shared<MockIntegratorCore>();

    //double on each call
    EXPECT_CALL(*mock_core, step(_, _, _, Eq(1), _)).Times(1).WillOnce(DoStepAndSetStepsize(2.));
    EXPECT_CALL(*mock_core, step(_, _, _, Eq(2), _)).Times(1).WillOnce(DoStepAndSetStepsize(4.));
    EXPECT_CALL(*mock_core, step(_, _, _, Eq(4), _)).Times(1).WillOnce(DoStepAndSetStepsize(8.));
    //last step of each advance call is smaller
    EXPECT_CALL(*mock_core, step(_, _, _, Eq(3), _)).Times(2).WillRepeatedly(DoStepAndSetStepsize(6.));
    //continue on the second advance call with updated step size, but last step of first advance is ignored
    EXPECT_CALL(*mock_core, step(_, _, _, Eq(8), _)).Times(1).WillOnce(DoStepAndSetStepsize(16.));

    auto f          = [](auto&& y, auto&& t, auto&& dydt) {};
    auto integrator = epi::OdeIntegrator(f, 0, Eigen::VectorXd::Constant(1, 0), 1.0, mock_core);
    integrator.advance(10.0);
    integrator.advance(21.0);
}
