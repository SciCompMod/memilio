#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <epidemiology/euler.h>

void init_vectors(std::vector<Eigen::VectorXd>& y, std::vector<Eigen::VectorXd>& sol, size_t n)
{
    y   = std::vector<Eigen::VectorXd>(n, Eigen::VectorXd::Constant(1, 0));
    sol = std::vector<Eigen::VectorXd>(n, Eigen::VectorXd::Constant(1, 0));
}

// Test for y'(t) = cos(t)
void integration_test(std::vector<Eigen::VectorXd>& y, std::vector<Eigen::VectorXd>& sol, size_t n,
                      double dt, double& err)
{

    sol[0][0]     = std::sin(0);
    sol[n - 1][0] = std::sin((n - 1) * dt);
    auto f = [](Eigen::VectorXd const& y, const double t, Eigen::VectorXd& dydt) { 
        dydt[0] = std::cos(t); 
    };

    double t = 0.;
    for (size_t i = 0; i < n - 1; i++) {
        sol[i + 1][0] = std::sin((i + 1) * dt);

        epi::EulerIntegratorCore().step(f, y[i], t, dt, y[i + 1]);

        printf("\n %.8f\t %.8f", y[i + 1][0], sol[i + 1][0]);
        // printf("\n approx: %.4e, sol: %.4e, error %.4e", y[i+1][0], sol[i+1][0], err);

        err += std::pow(std::abs(y[i + 1][0] - sol[i + 1][0]), 2.0);
    }
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

    printf("\n .%.8f. \n", dt);

    init_vectors(y, sol, n);
    integration_test(y, sol, n, dt, err);

    err = std::sqrt(err) / n;

    printf("\nFor n=%d the error is %.4e\n", (int)n, err);
}
