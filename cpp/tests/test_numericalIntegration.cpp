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

    Eigen::VectorXd f = Eigen::VectorXd::Constant(1, 0);

    sol[0][0]     = std::sin(0);
    sol[n - 1][0] = std::sin((n - 1) * dt);

    epi::EulerIntegrator euler([](Eigen::VectorXd const& y, const double t, Eigen::VectorXd& dydt) {
        dydt[0] = std::cos(t);
    });

    for (size_t i = 0; i < n - 1; i++) {
        sol[i + 1][0] = std::sin((i + 1) * dt);

        euler.step(y[i], t, dt, y[i + 1]);

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

    epi::RKIntegrator rkf45(sin_deriv, 1e-3, 1.0);
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

        rkf45.step(y[i], t_eval, dt, y[i + 1]); //

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

class FakeNonAdaptiveIntegrator : public epi::IntegratorBase
{
    public:
    FakeNonAdaptiveIntegrator(epi::DerivFunction f) : epi::IntegratorBase(f)
    {}
    virtual bool step(const Eigen::VectorXd& yt, double& t, double& dt, Eigen::VectorXd& ytp1) const override
    {
        t += dt;
        ytp1 = yt;
        return true;
    }
};

template <class FakeIntegrator>
class MockIntegrator : public epi::IntegratorBase
{
public:
    MockIntegrator(epi::DerivFunction f)
        : epi::IntegratorBase(f)
        , m_fake(f)
    {
        ON_CALL(*this, step).WillByDefault([this](auto&& yt, auto&& t, auto&& dt, auto&& ytp1) {
            return m_fake.step(yt, t, dt, ytp1);
        });
    }

    MOCK_METHOD(bool, step, (const Eigen::VectorXd& yt, double& t, double& dt, Eigen::VectorXd& ytp1), (const));

private:
    FakeIntegrator m_fake;
};

TEST(TestOdeIntegrate, integratorDoesTheRightNumberOfSteps)
{
    using testing::_;
    MockIntegrator<FakeNonAdaptiveIntegrator> mockIntegrator([](const auto& y, auto t, auto& dydt) {});
    EXPECT_CALL(mockIntegrator, step(_, _, _, _)).Times(100);
    std::vector<Eigen::VectorXd> result_x(1, Eigen::VectorXd::Constant(1, 0.0));
    auto result_t = epi::ode_integrate(0, 1, 1e-2, mockIntegrator, result_x);
    EXPECT_EQ(result_t.size(), 101);
    EXPECT_EQ(result_x.size(), 101);
}

TEST(TestOdeIntegrate, integratorStopsAtTMax)
{
    std::vector<Eigen::VectorXd> result(1, Eigen::VectorXd::Constant(1, 0));
    auto t =
        epi::ode_integrate(0, 2.34, 0.137, FakeNonAdaptiveIntegrator([](const auto& y, auto t, auto& dydt) {}), result);
    EXPECT_FLOAT_EQ(t.back(), 2.34);
}
