#include <epidemiology/utils/eigen_util.h>
#include <epidemiology/utils/stl_util.h>
#include "epidemiology/utils/compiler_diagnostics.h"
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
    IF_CONSTEXPR(TypeParam::IsRowMajor) {
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

TEST(TestEigenUtil, max)
{
    auto A = (Eigen::MatrixXi(2, 3) << -1, 2, -3, -4, 5, 6).finished();
    auto M = epi::max(A, Eigen::MatrixXi::Zero(2, 3));
    EXPECT_EQ(print_wrap(M), print_wrap((Eigen::MatrixXi(2, 3) << 0, 2, 0, 0, 5, 6).finished()));
}

TEST(TestRowMajorIterator, in_memory)
{
    auto m = Eigen::MatrixXd(3, 2);
    m << 0, 1, 2, 3, 4, 5;

    EXPECT_THAT(epi::make_range(epi::begin(m), epi::end(m)), testing::ElementsAre(0.0, 1.0, 2., 3., 4., 5.));
    EXPECT_THAT(epi::make_range(epi::cbegin(m), epi::cend(m)), testing::ElementsAre(0.0, 1.0, 2., 3., 4., 5.));
}

TEST(TestRowMajorIterator, expression)
{
    auto m = Eigen::MatrixXd::Constant(3, 2, 0.5);

    EXPECT_THAT(epi::make_range(epi::begin(m), epi::end(m)), testing::ElementsAre(0.5, 0.5, 0.5, 0.5, 0.5, 0.5));
    EXPECT_THAT(epi::make_range(epi::cbegin(m), epi::cend(m)), testing::ElementsAre(0.5, 0.5, 0.5, 0.5, 0.5, 0.5));
}

TEST(TestRowMajorIterator, operators)
{
    auto m = Eigen::MatrixXd(3, 2);
    m << 0, 1, 2, 3, 4, 5;

    auto b = epi::begin(m);
    auto e = epi::end(m);

    {
        auto it = b + 3;
        EXPECT_EQ(it, --(++it));
    }

    EXPECT_EQ(*(b + 2), 2.0);
    EXPECT_EQ(*(3 + b), 3.0);
    EXPECT_EQ(*(e - 2), 4.0);
    EXPECT_EQ(e - b, 6);

    EXPECT_TRUE(b < e);
    EXPECT_TRUE(b <= e);
    EXPECT_TRUE(b + 6 <= e);
    EXPECT_TRUE(e > b);
    EXPECT_TRUE(e >= b);
    EXPECT_TRUE(e - 6 >= b);
    
    EXPECT_EQ(b + 6, e);
    EXPECT_NE(b + 5, e);
}