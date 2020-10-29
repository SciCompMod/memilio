#include <epidemiology_io/twitter_migration_io.h>
#include <epidemiology/utils/logging.h>
#include "epidemiology/utils/eigen.h"

#include <gtest/gtest.h>

TEST(TestTwitter, compareWithSimpleFile)
{
    Eigen::MatrixXi test_matrix(4, 4);
    for (int i = 0; i < 4; i++) {
        test_matrix(i, i) = 0;
        for (int j = i + 1; j < 4; j++) {
            test_matrix(i, j) = i * 4 + j;
            test_matrix(j, i) = (i * 4 + j) * 10;
        }
    }

    std::fstream file;

    file.open("test_twitter.txt", std::ios::out);

    if (!file) {
        epi::log_error("File was not created");
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

    auto matrix_read = epi::read_migration("test_twitter.txt");

    ASSERT_EQ(test_matrix.rows(), matrix_read.rows());
    ASSERT_EQ(test_matrix.cols(), matrix_read.cols());
    for (int i = 0; i < test_matrix.rows(); i++) {
        for (int j = 0; j < test_matrix.cols(); j++) {
            ASSERT_EQ(test_matrix(i, j), matrix_read(i, j));
        }
    }
}
