/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Wadim Koslow
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
#include "memilio/io/mobility_io.h"
#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"
#include "matchers.h"

#include <gtest/gtest.h>

TEST(TestReadMobility, readFormatted)
{
    Eigen::MatrixXd test_matrix(4, 4);
    for (int i = 0; i < 4; i++) {
        test_matrix(i, i) = 0;
        for (int j = i + 1; j < 4; j++) {
            test_matrix(i, j) = i * 4 + j;
            test_matrix(j, i) = (i * 4 + j) * 10;
        }
    }

    std::fstream file;

    file.open("test_mobility.txt", std::ios::out);

    if (!file) {
        mio::log_error("File was not created");
    }
    else {
        file << "from_str\tto_str\tfrom_rs\tto_rs\tcount_abs\n";
        file << "Narnia\tHogwarts\t150\t42\t" + std::to_string(test_matrix(3, 2)) + "\n";
        file << "Middle-Earth\tWesteros\t12\t7\t" + std::to_string(test_matrix(1, 0)) + "\n";
        file << "Narnia\tMiddle-Earth\t150\t12\t" + std::to_string(test_matrix(3, 1)) + "\n";
        file << "Narnia\tWesteros\t150\t7\t" + std::to_string(test_matrix(3, 0)) + "\n";
        file << "Hogwarts\tNarnia\t42\t150\t" + std::to_string(test_matrix(2, 3)) + "\n";
        file << "Hogwarts\tMiddle-Earth\t42\t12\t" + std::to_string(test_matrix(2, 1)) + "\n";
        file << "Hogwarts\tWesteros\t42\t7\t" + std::to_string(test_matrix(2, 0)) + "\n";
        file << "MiddleEarth\tHogwarts\t12\t42\t" + std::to_string(test_matrix(1, 2)) + "\n";
        file << "MiddleEarth\tNarnia\t12\t150\t" + std::to_string(test_matrix(1, 3)) + "\n";
        file << "Westeros\tNarnia\t7\t150\t" + std::to_string(test_matrix(0, 3)) + "\n";
        file << "Westeros\tHogwarts\t7\t42\t" + std::to_string(test_matrix(0, 2)) + "\n";
        file << "Westeros\tMiddle-Earth\t7\t12\t" + std::to_string(test_matrix(0, 1)) + "\n";

        file.close();
    }

    auto matrix_read = mio::read_mobility_formatted("test_mobility.txt");
    ASSERT_TRUE(matrix_read);
    ASSERT_EQ(test_matrix.rows(), matrix_read.value().rows());
    ASSERT_EQ(test_matrix.cols(), matrix_read.value().cols());
    ASSERT_EQ(print_wrap(test_matrix), print_wrap(matrix_read.value()));
}

TEST(TestReadMobility, readPlain)
{
    Eigen::MatrixXd test_matrix(6, 6);

    test_matrix(0, 0) = 0.4413;
    test_matrix(0, 1) = 0.4504;
    test_matrix(0, 2) = 1.2383;
    test_matrix(0, 3) = 0.8033;
    test_matrix(0, 4) = 0.0494;
    test_matrix(0, 5) = 0.0017;
    test_matrix(1, 0) = 0.0485;
    test_matrix(1, 1) = 0.7616;
    test_matrix(1, 2) = 0.6532;
    test_matrix(1, 3) = 1.1614;
    test_matrix(1, 4) = 0.0256;
    test_matrix(1, 5) = 0.0013;
    test_matrix(2, 0) = 0.1800;
    test_matrix(2, 1) = 0.1795;
    test_matrix(2, 2) = 0.8806;
    test_matrix(2, 3) = 0.6413;
    test_matrix(2, 4) = 0.0429;
    test_matrix(2, 5) = 0.0032;
    test_matrix(3, 0) = 0.0495;
    test_matrix(3, 1) = 0.2639;
    test_matrix(3, 2) = 0.5189;
    test_matrix(3, 3) = 0.8277;
    test_matrix(3, 4) = 0.0679;
    test_matrix(3, 5) = 0.0014;
    test_matrix(4, 0) = 0.0087;
    test_matrix(4, 1) = 0.0394;
    test_matrix(4, 2) = 0.1417;
    test_matrix(4, 3) = 0.3834;
    test_matrix(4, 4) = 0.7064;
    test_matrix(4, 5) = 0.0447;
    test_matrix(5, 0) = 0.0292;
    test_matrix(5, 1) = 0.0648;
    test_matrix(5, 2) = 0.1248;
    test_matrix(5, 3) = 0.4179;
    test_matrix(5, 4) = 0.3497;
    test_matrix(5, 5) = 0.1544;

    auto matrix_read = mio::read_mobility_plain(get_test_data_file_path("contacts.txt"));

    ASSERT_TRUE(matrix_read);
    ASSERT_EQ(test_matrix.rows(), matrix_read.value().rows());
    ASSERT_EQ(test_matrix.cols(), matrix_read.value().cols());
    ASSERT_EQ(print_wrap(test_matrix), print_wrap(matrix_read.value()));
}
