/*
* Copyright (C) 2020-2026 MEmilio
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
#include "memilio/epidemiology/contact_matrix.h"
#include "matchers.h"
#include "gtest/gtest.h"

TEST(TestContactMatrix, initZero)
{
    mio::ContactMatrix<double> cm(Eigen::Index(3));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(mio::SimulationTime<double>(-1e5))), print_wrap(Eigen::MatrixXd::Zero(3, 3)));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(mio::SimulationTime<double>(0))), print_wrap(Eigen::MatrixXd::Zero(3, 3)));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(mio::SimulationTime<double>(1e-32))),
              print_wrap(Eigen::MatrixXd::Zero(3, 3)));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(mio::SimulationTime<double>(1e5))), print_wrap(Eigen::MatrixXd::Zero(3, 3)));
}

TEST(TestContactMatrix, initBaseAndMin)
{
    auto B = (Eigen::MatrixXd(2, 2) << 1, 2, 3, 4).finished();
    auto M = Eigen::MatrixXd::Constant(2, 2, 0.1);
    mio::ContactMatrix<double> cm(B, M);
    EXPECT_EQ(print_wrap(cm.get_matrix_at(mio::SimulationTime<double>(-1e5))), print_wrap(B));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(mio::SimulationTime<double>(0))), print_wrap(B));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(mio::SimulationTime<double>(1e-32))), print_wrap(B));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(mio::SimulationTime<double>(1e5))), print_wrap(B));
}

TEST(TestContactMatrix, dampings)
{
    auto B = (Eigen::MatrixXd(2, 2) << 1, 2, 3, 4).finished();
    auto M = Eigen::MatrixXd::Constant(2, 2, 0.1);
    mio::ContactMatrix<double> cm(B, M);
    auto D1 = 0.25;
    auto D2 = (Eigen::MatrixXd(2, 2) << 0.0, 0.75, 0.5, 0.25).finished();
    cm.add_damping(D1, mio::DampingLevel(7), mio::DampingType(3), mio::SimulationTime<double>(0.5));
    cm.add_damping(D2, mio::DampingLevel(7), mio::DampingType(2), mio::SimulationTime<double>(2.0));

    EXPECT_EQ(cm.get_dampings().size(), 2);
    EXPECT_THAT(cm.get_dampings(),
                testing::ElementsAre(mio::SquareDamping<double>(D1, mio::DampingLevel(7), mio::DampingType(3),
                                                                mio::SimulationTime<double>(0.5), Eigen::Index(2)),
                                     mio::SquareDamping<double>(D2, mio::DampingLevel(7), mio::DampingType(2),
                                                                mio::SimulationTime<double>(2.0))));

    EXPECT_EQ(print_wrap(cm.get_matrix_at(mio::SimulationTime<double>(-1e5))), print_wrap(B));
    EXPECT_EQ(print_wrap(cm.get_matrix_at(mio::SimulationTime<double>(-0.5))), print_wrap(B));
    EXPECT_THAT(print_wrap(cm.get_matrix_at(mio::SimulationTime<double>(0.5 + 1e-32))), MatrixNear(B - D1 * (B - M)));
    EXPECT_THAT(print_wrap(cm.get_matrix_at(mio::SimulationTime<double>(1e5))),
                MatrixNear(B - ((D1 + D2.array()) * (B - M).array()).matrix()));
}

TEST(TestContactMatrixGroup, sum)
{
    mio::ContactMatrixGroup<double> cmg(3, 2);
    cmg[0] = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(3, 3, 1.0));
    cmg[1] = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(3, 3, 2.0));
    cmg[2] = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(3, 3, 3.0));
    cmg.add_damping(0.5, mio::DampingLevel(3), mio::DampingType(1), mio::SimulationTime<double>(1.0));

    EXPECT_THAT(print_wrap(cmg.get_matrix_at(mio::SimulationTime<double>(0.0))),
                MatrixNear(Eigen::MatrixXd::Constant(3, 3, 6.0)));
    EXPECT_THAT(print_wrap(cmg.get_matrix_at(mio::SimulationTime<double>(1.0))),
                MatrixNear(Eigen::MatrixXd::Constant(3, 3, 3.0)));
}
