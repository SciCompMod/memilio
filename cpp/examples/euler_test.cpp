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
#include "memilio/math/euler.h"

#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <iostream>
#include <cmath>

void init_vectors(std::vector<Eigen::VectorXd>& y, std::vector<Eigen::VectorXd>& sol, size_t n)
{
    y   = std::vector<Eigen::VectorXd>(n, Eigen::VectorXd::Constant(1, 0));
    sol = std::vector<Eigen::VectorXd>(n, Eigen::VectorXd::Constant(1, 0));
}

// Test for y'(t) = cos(t)
void integration_test(std::vector<Eigen::VectorXd>& y, std::vector<Eigen::VectorXd>& sol, size_t n, double dt,
                      double& err)
{

    sol[0][0]     = std::sin(0);
    sol[n - 1][0] = std::sin((n - 1) * dt);
    auto f        = [](auto&& /*y*/, auto&& t, auto&& dydt) {
        dydt[0] = std::cos(t);
    };

    double t = 0.;
    for (size_t i = 0; i < n - 1; i++) {
        sol[i + 1][0] = std::sin((i + 1) * dt);

        mio::EulerIntegratorCore().step(f, y[i], t, dt, y[i + 1]);

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
