#include "epidemiology/utils/matrix_shape.h"
#include "gtest/gtest.h"

TEST(TestMatrixShape, Rect)
{
    EXPECT_EQ(epi::RectMatrixShape(3, 2).rows(), 3);
    EXPECT_EQ(epi::RectMatrixShape(3, 2).cols(), 2);
    EXPECT_EQ(epi::RectMatrixShape::get_shape_of(Eigen::MatrixXd::Zero(2, 3)), epi::RectMatrixShape(2, 3));
}

TEST(TestMatrixShape, Square)
{
    EXPECT_EQ(epi::SquareMatrixShape(3).rows(), 3);
    EXPECT_EQ(epi::SquareMatrixShape(3).cols(), 3);
    EXPECT_EQ(epi::SquareMatrixShape(3).size(), 3);
    EXPECT_EQ(epi::SquareMatrixShape::get_shape_of(Eigen::MatrixXd::Zero(3, 3)), epi::SquareMatrixShape(3));
}

TEST(TestMatrixShape, Vector)
{
    EXPECT_EQ(epi::ColumnVectorShape(3).rows(), 3);
    EXPECT_EQ(epi::ColumnVectorShape(3).cols(), 1);
    EXPECT_EQ(epi::ColumnVectorShape(3).size(), 3);
    EXPECT_EQ(epi::ColumnVectorShape::get_shape_of(Eigen::VectorXd::Zero(3)), epi::ColumnVectorShape(3));
}