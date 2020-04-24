#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <epidemiology/euler.h>


template <typename T>
void init_vectors(std::vector<std::vector<T> > &y, std::vector<std::vector<T> > &sol, std::vector<T> &f, size_t n)
{
    y = std::vector<std::vector<T> >(n, std::vector<T>(1, 0));
    sol = std::vector<std::vector<T> >(n, std::vector<T>(1, 0));

    f = std::vector<T>(1, 0);
}

// Test for y'(t) = cos(t)
template <typename T>
void integration_test(std::vector<std::vector<T> > &y, std::vector<std::vector<T> > &sol, std::vector<T> &f, size_t n, T dt, T& err)
{

    sol[0][0] = std::sin(0);
    sol[n - 1][0] = std::sin((n - 1) * dt);


    for(size_t i = 0; i < n - 1; i++) {
        sol[i + 1][0] = std::sin((i + 1) * dt);

        f[0] = std::cos(i * dt);

        explicit_euler(y[i], f, dt, y[i + 1]); //

printf("\n %.8f\t %.8f", y[i + 1][0], sol[i + 1][0]);
        // printf("\n approx: %.4e, sol: %.4e, error %.4e", y[i+1][0], sol[i+1][0], err);

        err += std::pow(std::abs(y[i + 1][0] - sol[i + 1][0]), 2.0);

    }

}

int main()
{
    std::vector<std::vector<double> > y;
    std::vector<std::vector<double> > sol;

    std::vector<double> f;

    const double pi = std::acos(-1);

    size_t n = 10;
    double t0 = 0;
    double tmax = 2 * pi;
    double dt = (tmax-t0)/n;
    double err = 0;

    printf("\n .%.8f. \n", dt);

    init_vectors(y, sol, f, n);
    integration_test(y, sol, f, n, dt, err);

    err = std::sqrt(err) / n;

    printf("\nFor n=%d the error is %.4e\n", (int)n, err);
}
