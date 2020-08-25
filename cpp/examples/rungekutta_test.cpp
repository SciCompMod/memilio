#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <epidemiology/euler.h>
#include <epidemiology/adapt_rk.h>

void init_vectors(std::vector<Eigen::VectorXd>& y, std::vector<Eigen::VectorXd>& sol, size_t n)
{
    y   = std::vector<Eigen::VectorXd>(n, Eigen::VectorXd::Constant(1, 0));
    sol = std::vector<Eigen::VectorXd>(n, Eigen::VectorXd::Constant(1, 0));
}

// Test for y'(t) = cos(t)
void integration_test(std::vector<Eigen::VectorXd>& y, std::vector<Eigen::VectorXd>& sol, size_t& n, double t, double dt,
                      const double tmax, double& err)
{
    auto sine_deriv = [](auto&& y, auto&& t, auto&& dydt) {
        dydt[0] = std::cos(t);
    };

    epi::RKIntegratorCore rkf45(1e-3, 1.0);
    rkf45.set_abs_tolerance(1e-7);
    rkf45.set_rel_tolerance(1e-7);

    sol[0][0] = std::sin(0);

    std::vector<double> f = std::vector<double>(1, 0);
    size_t i         = 0;
    double t_eval    = t;
    // printf("\n t: %.8f\t sol %.8f\t rkf %.8f", t, sol[0][0], y[0][0]);

    while (t_eval - tmax < 1e-10) {

        if (i + 1 >= sol.size()) {
            sol.push_back(Eigen::VectorXd::Constant(1, 0));
            y.push_back(Eigen::VectorXd::Constant(1, 0));
        }

        double dt_old = dt;

        rkf45.step(sine_deriv, y[i], t_eval, dt, y[i + 1]); //

        sol[i + 1][0] = std::sin(t_eval);

        // printf("\n t: %.8f (dt %.8f)\t sol %.8f\t rkf %.8f", t_eval, dt, sol[i + 1][0], y[i + 1][0]);
        // printf("\n approx: %.4e, sol: %.4e, error %.4e", y[i+1][0], sol[i+1][0], err);

        err += std::pow(std::abs(y[i + 1][0] - sol[i + 1][0]), 2.0);
        i++;
    }

    n = i;
}

int main()
{
    std::vector<Eigen::VectorXd> y;
    std::vector<Eigen::VectorXd> sol;

    const double pi = std::acos(-1);

    size_t n    = 10;
    double t0   = 0;
    double tmax = 2 * pi;
    double dt   = (tmax - t0) / n;
    double err  = 0;

    init_vectors(y, sol, n);

    integration_test(y, sol, n, t0, dt, tmax, err);

    err = std::sqrt(err) / n;

    printf("\nFor n=%d the error is %.4e\n", (int)n, err);
}
