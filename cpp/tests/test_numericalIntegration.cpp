#include <gtest/gtest.h>
#include <epidemiology/euler.h>
#include <epidemiology/adapt_rk.h>

#include <string>
#include <vector>
#include <fstream>
#include <ios>
#include <cmath>

void sin_deriv(std::vector<double> const &y, const double t, std::vector<double> &dydt)
{
    dydt[0] = std::cos(t);
}

class TestVerifyNumericalIntegrator : public testing::Test
{
protected:

    void SetUp() override
    {
        t = 0.;
        tmax = 2 * std::acos(-1); // 2PI
        err = 0;
    }

public:

    std::vector<std::vector<double> > y;
    std::vector<std::vector<double> > sol;
   
    double t;
    double tmax;
    size_t n;
    double dt;
    double err;

};


TEST_F(TestVerifyNumericalIntegrator, euler_sine)
{
    n = 1000;
    dt = (tmax-t)/n;
    y = std::vector<std::vector<double>>(n, std::vector<double>(1, 0));
    sol = std::vector<std::vector<double>>(n, std::vector<double>(1, 0));
    
    std::vector<double> f = std::vector<double>(1, 0);

    sol[0][0] = std::sin(0);
    sol[n - 1][0] = std::sin((n - 1) * dt);


    for(size_t i = 0; i < n - 1; i++) {
        sol[i + 1][0] = std::sin((i + 1) * dt);

        f[0] = std::cos(i * dt);

        explicit_euler(y[i], f, dt, y[i + 1]); //

        // printf("\n %.8f\t %.8f ", y[i + 1][0], sol[i + 1][0]);

        err += std::pow(std::abs(y[i + 1][0] - sol[i + 1][0]), 2.0);

    }

    err = std::sqrt(err) / n;

    EXPECT_NEAR(err, 0.0, 1e-3);

}



TEST_F(TestVerifyNumericalIntegrator, runge_kutta_fehlberg45_sine)
{

    n = 10;
    dt = (tmax-t)/n;
    y = std::vector<std::vector<double>>(n, std::vector<double>(1, 0));
    sol = std::vector<std::vector<double>>(n, std::vector<double>(1, 0));

    RKIntegrator<double> rkf45(sin_deriv, 1e-3, 1.0);
    rkf45.set_abs_tolerance(1e-7);
    rkf45.set_rel_tolerance(1e-7);

    sol[0][0] = std::sin(0);

    size_t i = 0; 
    double t_eval = t;
    // printf("\n t: %.8f\t sol %.8f\t rkf %.8f", t, sol[0][0], y[0][0]);
    while(t_eval-tmax < 1e-10) {

        if(i+1 >= sol.size())
        {
            sol.push_back(std::vector<double>(1, 0));
            y.push_back(std::vector<double>(1, 0));
        }

        double dt_old = dt;

        rkf45.step(y[i], t_eval, dt, y[i + 1]); //

        sol[i + 1][0] = std::sin(t_eval);

        // printf("\n t: %.8f (dt %.8f)\t sol %.8f\t rkf %.8f", t_eval, dt, sol[i + 1][0], y[i + 1][0]);
        // printf("\n approx: %.4e, sol: %.4e, error %.4e", y[i+1][0], sol[i+1][0], err);

        err += std::pow(std::abs(y[i + 1][0] - sol[i + 1][0]), 2.0);
        i++;
    }

    n=i;

    err = std::sqrt(err) / n;

    EXPECT_NEAR(err, 0.0, 1e-7);

}
