#include <epidemiology/utils/eigen_util.h>
#include "matchers.h"
#include <gtest/gtest.h>

namespace Eigen
{
void PrintTo(const Eigen::VectorXd& m, std::ostream* os)
{
    (*os) << m;
}
} // namespace Eigen

template <class M>
struct TestEigenUtilMatrix : public ::testing::Test {
};
using MatrixTypes =
    ::testing::Types<Eigen::MatrixXd, Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>;
TYPED_TEST_SUITE(TestEigenUtilMatrix, MatrixTypes);
TYPED_TEST(TestEigenUtilMatrix, slice)
{
    TypeParam A(5, 4);
    A << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19;

    TypeParam B(3, 2);
    B << 4, 6, 8, 10, 12, 14;

    ASSERT_EQ(print_wrap(epi::slice(A, {1, 3, 1}, {0, 2, 2})), print_wrap(B));
}

template <class V>
struct TestEigenUtilVector : public ::testing::Test {
};
using VectorTypes = ::testing::Types<Eigen::VectorXd, Eigen::RowVectorXd>;
TYPED_TEST_SUITE(TestEigenUtilVector, VectorTypes);
TYPED_TEST(TestEigenUtilVector, slice)
{
    TypeParam A(12);
    A << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11;

    TypeParam B(4);
    B << 1, 4, 7, 10;

    ASSERT_EQ(print_wrap(epi::slice(A, {1, 4, 3})), print_wrap(B));
}

TYPED_TEST(TestEigenUtilMatrix, reshape)
{
    TypeParam A(2, 3), B(1, 6), C(3, 2), D(6, 1);
    A << 0, 1, 2, 3, 4, 5;
    if (TypeParam::IsRowMajor) {
        B << 0, 1, 2, 3, 4, 5;
        C << 0, 1, 2, 3, 4, 5;
        D << 0, 1, 2, 3, 4, 5;
    }
    else {
        B << 0, 3, 1, 4, 2, 5;
        C << 0, 4, 3, 2, 1, 5;
        D << 0, 3, 1, 4, 2, 5;
    }

    EXPECT_EQ(print_wrap(epi::reshape(A, 1, 6)), print_wrap(B));
    EXPECT_EQ(print_wrap(epi::reshape(A, 3, 2)), print_wrap(C));
}