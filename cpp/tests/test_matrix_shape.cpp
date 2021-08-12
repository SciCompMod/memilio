/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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