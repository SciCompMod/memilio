/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler
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
#include "load_test_data.h"
#include "ide_sir/infection_state.h"
#include "ide_sir/model.h"
#include "ide_sir/parameters.h"
#include "ide_sir/simulation.h"
#include "memilio/config.h"
#include <gtest/gtest.h>
#include <vector>

// Check that the Gregory weights are set correctly for orders 1, 2 and 3.
TEST(IdeSir, checkGregoryWeights)
{

    // Define t_max. Below, we will test the Gregory weights for n=gregory_order, ..., n_max.
    size_t n_max = 6;

    // Define Gregory orders that we test.
    std::vector<size_t> gregory_orders = {2, 3};

    for (size_t gregory_order : gregory_orders) {

        // Here, we define the matrices for Sigma and the vector corresponding to Omega explicitly for n up to n_max. We
        // will then compare this to our implementation in get_gregoryweights(), where we implemented the matrices and
        // vectors in a reduced form independent of n_max.
        std::vector<std::vector<ScalarType>> gregoryWeights_sigma_expected;
        std::vector<ScalarType> gregoryWeights_omega_expected;
        std::cout << "gregory order: " << gregory_order << std::endl;
        switch (gregory_order) {
        case 1:
            // TODO
            break;
        case 2:
            std::cout << "case 2 \n";
            gregoryWeights_sigma_expected = {{5. / 12., 14. / 12.},
                                             {5. / 12., 13. / 12.},
                                             {5. / 12., 13. / 12.},
                                             {5. / 12., 13. / 12.},
                                             {5. / 12., 13. / 12.}};
            gregoryWeights_omega_expected = {5. / 12., 13. / 12., 12. / 12., 12. / 12., 12. / 12.};
            break;
        case 3:
            std::cout << "case 3 \n";
            gregoryWeights_sigma_expected = {{9. / 24., 27. / 24., 27. / 24.},
                                             {9. / 24., 28. / 24., 22. / 24.},
                                             {9. / 24., 28. / 24., 23. / 24.},
                                             {9. / 24., 28. / 24., 23. / 24.}};
            gregoryWeights_omega_expected = {9. / 24., 28. / 24., 23. / 24., 24. / 24.};
            break;
        }
        // Get Gregory weights from implementation in get_gregoryweights().
        std::vector<Eigen::MatrixX<ScalarType>> vec_gregoryweights = mio::isir::get_gregoryweights(gregory_order);
        Eigen::MatrixX<ScalarType> gregoryWeights_sigma            = vec_gregoryweights[0];
        Eigen::MatrixX<ScalarType> gregoryWeights_omega            = vec_gregoryweights[1];

        size_t row_index, column_index;
        for (size_t n = gregory_order; n <= n_max; n++) {
            std::cout << "n: " << n << std::endl;
            for (size_t j = 0; j < gregory_order; j++) {
                std::cout << "j: " << j << std::endl;

                // Check sigma.

                // Compute row_index and column_index according to n and j as in the implementation of sum_part1.
                // If n <= gregory_order + gregory_order - 2, then the row_index corresponds to n - gregory_order.
                if (n <= gregory_order + gregory_order - 2) {
                    row_index = n - gregory_order;
                }
                // Else, for n >= m_gregory_order - 1, the entries in gregoryWeights_sigma do not change anymore and the
                // corresponding row_index is given by gregory_order - 1.
                else {
                    row_index = gregory_order - 1;
                }
                // The column index only depends on the current index of the sum j.
                column_index = j;

                // std::cout << "lhs: " << gregoryWeights_sigma(row_index, column_index) << std::endl;
                // std::cout << "rhs: " << gregoryWeights_sigma_expected[n - gregory_order][j] << std::endl;
                EXPECT_EQ(gregoryWeights_sigma(row_index, column_index),
                          gregoryWeights_sigma_expected[n - gregory_order][j]);
            }
        }

        size_t weight_index;
        for (size_t n = gregory_order; n <= n_max; n++) {
            // std::cout << "n: " << n << std::endl;
            for (size_t j = gregory_order; j <= n; j++) {
                // std::cout << "j: " << j << std::endl;
                // Check omega.

                // Compute weight_index according to n and j as in the implementation of sum_part2.
                // Depending on the Gregory order and the current time step, we determine the required weight index of gregoryWeights_omega.
                // This is necessary because we implemented gregoryWeights_omega in a reduced way.
                if (n - j <= gregory_order) {
                    weight_index = n - j;
                }
                else {
                    weight_index = gregory_order;
                }

                mio::unused(weight_index);

                // std::cout << "n-j: " << n - j << std::endl;

                EXPECT_EQ(gregoryWeights_omega(weight_index), gregoryWeights_omega_expected[n - j]);
            }

            // Define weight_indices that we want to check the values for in gregoryWeights_omega.
            std::vector<size_t> weight_indices = {0, 1, 2, 3};
        }
    }
}
