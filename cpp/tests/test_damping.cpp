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
#include "memilio/epidemiology/damping.h"
#include "matchers.h"
#include <gtest/gtest.h>

TEST(TestDampings, initZero)
{
    mio::Dampings<double, mio::Damping<double, mio::RectMatrixShape<double>>> dampings(3, 2);
    EXPECT_EQ(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(-1e5))),
              print_wrap(Eigen::MatrixXd::Zero(3, 2)));
    EXPECT_EQ(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(0))),
              print_wrap(Eigen::MatrixXd::Zero(3, 2)));
    EXPECT_EQ(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(1e-32))),
              print_wrap(Eigen::MatrixXd::Zero(3, 2)));
    EXPECT_EQ(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(1e5))),
              print_wrap(Eigen::MatrixXd::Zero(3, 2)));
}

TEST(TestDampings, dampingsOnDifferentLevels)
{
    mio::Dampings<double, mio::Damping<double, mio::RectMatrixShape<double>>> dampings(2, 2);
    auto D1 = 0.25;
    auto D2 = (Eigen::MatrixXd(2, 2) << 0.25, 0.5, 0.75, 1).finished();
    dampings.add(D1, mio::DampingLevel(7), mio::DampingType(3), mio::SimulationTime<double>(0.5));
    dampings.add(D2, mio::DampingLevel(13), mio::DampingType(3), mio::SimulationTime<double>(2.0));
    EXPECT_EQ(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(-1e5))),
              print_wrap(Eigen::MatrixXd::Zero(2, 2)));
    EXPECT_EQ(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(-0.5))),
              print_wrap(Eigen::MatrixXd::Zero(2, 2)));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(0.5 + 1e-32))),
                MatrixNear(Eigen::MatrixXd::Constant(2, 2, D1)));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(1e5))),
                MatrixNear((D1 + D2.array() - D1 * D2.array()).matrix()));
}

TEST(TestDampings, dampingsOnSameLevel)
{
    mio::Dampings<double, mio::Damping<double, mio::SquareMatrixShape<double>>> dampings(2);
    auto D1 = 0.25;
    auto D2 = (Eigen::MatrixXd(2, 2) << 0.0, 0.25, 0.5, 0.75).finished();
    dampings.add(D1, mio::DampingLevel(-2), mio::DampingType(0), mio::SimulationTime<double>(0.5));
    dampings.add(D2, mio::DampingLevel(-2), mio::DampingType(1), mio::SimulationTime<double>(2.0));
    EXPECT_EQ(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(-1e5))),
              print_wrap(Eigen::MatrixXd::Zero(2, 2)));
    EXPECT_EQ(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(-0.5))),
              print_wrap(Eigen::MatrixXd::Zero(2, 2)));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(0.5 + 1e-32))),
                MatrixNear(Eigen::MatrixXd::Constant(2, 2, D1)));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(1e5))),
                MatrixNear((D1 + D2.array()).matrix()));
}

TEST(TestDampings, dampingsAtTheSameTime)
{
    mio::Dampings<double, mio::Damping<double, mio::SquareMatrixShape<double>>> dampings(2);
    auto D1 = 0.25;
    auto D2 = (Eigen::MatrixXd(2, 2) << 0.0, 0.25, 0.5, 0.75).finished();
    dampings.add(D1, mio::DampingLevel(-2), mio::DampingType(0), mio::SimulationTime<double>(0.5));
    dampings.add(D2, mio::DampingLevel(-2), mio::DampingType(1), mio::SimulationTime<double>(0.5));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(-0.5))),
                MatrixNear(Eigen::MatrixXd::Zero(2, 2)));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(0.5 + 1e-5))),
                MatrixNear((D1 + D2.array()).matrix()));
}

TEST(TestDampings, dampingOfSameType)
{
    mio::Dampings<double, mio::Damping<double, mio::SquareMatrixShape<double>>> dampings(2);
    auto D1 = 0.25;
    auto D2 = (Eigen::MatrixXd(2, 2) << 0.0, 0.25, 0.5, 0.75).finished();
    dampings.add(D1, mio::DampingLevel(123), mio::DampingType(5), mio::SimulationTime<double>(0.5));
    dampings.add(D2, mio::DampingLevel(123), mio::DampingType(5), mio::SimulationTime<double>(2.0));
    EXPECT_EQ(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(-1e5))),
              print_wrap(Eigen::MatrixXd::Zero(2, 2)));
    EXPECT_EQ(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(-0.5))),
              print_wrap(Eigen::MatrixXd::Zero(2, 2)));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(0.5 + 1e-32))),
                MatrixNear(Eigen::MatrixXd::Constant(2, 2, D1)));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(1e5))), MatrixNear(D2));
}

TEST(TestDampings, dampingsNegative)
{
    mio::Dampings<double, mio::Damping<double, mio::RectMatrixShape<double>>> dampings(2, 2);
    auto D1 = -4.25;
    auto D2 = (Eigen::MatrixXd(2, 2) << 0.25, -0.5, 0.75, -1).finished();
    dampings.add(D1, mio::DampingLevel(7), mio::DampingType(3), mio::SimulationTime<double>(0.5));
    dampings.add(D2, mio::DampingLevel(13), mio::DampingType(3), mio::SimulationTime<double>(2.0));
    EXPECT_EQ(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(-1e5))),
              print_wrap(Eigen::MatrixXd::Zero(2, 2)));
    EXPECT_EQ(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(-0.5))),
              print_wrap(Eigen::MatrixXd::Zero(2, 2)));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(0.5 + 1e-32))),
                MatrixNear(Eigen::MatrixXd::Constant(2, 2, D1)));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(1e5))),
                MatrixNear((D1 + D2.array() - D1 * D2.array()).matrix()));
}

TEST(TestDampings, dampingsCombined)
{
    mio::Dampings<double, mio::Damping<double, mio::SquareMatrixShape<double>>> dampings(2);
    auto D1 = 0.25;
    auto D2 = (Eigen::MatrixXd(2, 2) << 0.1, 0.1, 0.1, 0.1).finished();
    auto D3 = (Eigen::MatrixXd(2, 2) << 0.0, 0.25, 0.5, 0.75).finished();
    auto D4 = 0.5;
    //add dampings out of order to check sorting
    dampings.add(D2, mio::DampingLevel(7), mio::DampingType(2), mio::SimulationTime<double>(0.0));
    dampings.add(D1, mio::DampingLevel(123), mio::DampingType(5), mio::SimulationTime<double>(-2.0));
    dampings.add(D4, mio::DampingLevel(123), mio::DampingType(5), mio::SimulationTime<double>(3.0));
    dampings.add(D3, mio::DampingLevel(7), mio::DampingType(3), mio::SimulationTime<double>(1.5));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(-1e5))),
                print_wrap(Eigen::MatrixXd::Zero(2, 2)));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(-1.0))),
                MatrixNear(Eigen::MatrixXd::Constant(2, 2, D1)));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(0.2))),
                MatrixNear((D1 + D2.array() - D1 * D2.array()).matrix()));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(2.0))),
                MatrixNear((D1 + D2.array() + D3.array() - D1 * (D2 + D3).array()).matrix()));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(1e45))),
                MatrixNear((D4 + D2.array() + D3.array() - D4 * (D2 + D3).array()).matrix()));
}

TEST(TestDampings, smoothTransitions)
{
    mio::Dampings<double, mio::Damping<double, mio::ColumnVectorShape<double>>> dampings(2);
    auto D1 = 0.25;
    auto D2 = (Eigen::VectorXd(2) << 0.1, 0.1).finished();
    dampings.add(D1, mio::DampingLevel(123), mio::DampingType(5), mio::SimulationTime<double>(-2.0));
    dampings.add(D2, mio::DampingLevel(1), mio::DampingType(10), mio::SimulationTime<double>(1.5));

    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(-2.5))),
                MatrixNear((dampings.get_matrix_at(mio::SimulationTime<double>(-3.)) +
                            dampings.get_matrix_at(mio::SimulationTime<double>(-2.))) /
                           2));
    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(1.0))),
                MatrixNear((dampings.get_matrix_at(mio::SimulationTime<double>(0.5)) +
                            dampings.get_matrix_at(mio::SimulationTime<double>(1.5))) /
                           2));
}

TEST(TestDampings, automatic_cache_update)
{
    mio::Dampings<double, mio::Damping<double, mio::ColumnVectorShape<double>>> dampings(2);
    auto D1 = 0.25;
    dampings.set_automatic_cache_update(false);
    dampings.add(D1, mio::DampingLevel(1), mio::DampingType(2), mio::SimulationTime<double>(1.0));

#ifndef NDEBUG
    EXPECT_DEATH(dampings.get_matrix_at(mio::SimulationTime<double>(2.0)),
                 "Cache is not current\\. Did you disable the automatic cache update\\?");
#endif

    dampings.set_automatic_cache_update(true);

    EXPECT_THAT(print_wrap(dampings.get_matrix_at(mio::SimulationTime<double>(2.0))),
                MatrixNear((Eigen::VectorXd(2) << 0.25, 0.25).finished()));
}
