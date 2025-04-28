/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele
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
#include "memilio/math/matrix_shape.h"
#include "gtest/gtest.h"

TEST(TestMatrixShape, Rect)
{
    EXPECT_EQ(mio::RectMatrixShape(3, 2).rows(), 3);
    EXPECT_EQ(mio::RectMatrixShape(3, 2).cols(), 2);
    EXPECT_EQ(mio::RectMatrixShape::get_shape_of(Eigen::MatrixXd::Zero(2, 3)), mio::RectMatrixShape(2, 3));
}

TEST(TestMatrixShape, Square)
{
    EXPECT_EQ(mio::SquareMatrixShape(3).rows(), 3);
    EXPECT_EQ(mio::SquareMatrixShape(3).cols(), 3);
    EXPECT_EQ(mio::SquareMatrixShape(3).size(), 3);
    EXPECT_EQ(mio::SquareMatrixShape::get_shape_of(Eigen::MatrixXd::Zero(3, 3)), mio::SquareMatrixShape(3));
}

TEST(TestMatrixShape, Vector)
{
    EXPECT_EQ(mio::ColumnVectorShape(3).rows(), 3);
    EXPECT_EQ(mio::ColumnVectorShape(3).cols(), 1);
    EXPECT_EQ(mio::ColumnVectorShape(3).size(), 3);
    EXPECT_EQ(mio::ColumnVectorShape::get_shape_of(Eigen::VectorXd::Zero(3)), mio::ColumnVectorShape(3));
}
