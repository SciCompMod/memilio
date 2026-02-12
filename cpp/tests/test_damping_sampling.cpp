/*
* Copyright (C) 2020-2026 MEmilio
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
#include "memilio/epidemiology/damping_sampling.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "matchers.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>

TEST(TestDampingSampling, apply)
{
    auto ds =
        std::vector<mio::DampingSampling<double>>{mio::DampingSampling<double>{mio::UncertainValue<double>(0.5),
                                                                               mio::DampingLevel(0),
                                                                               mio::DampingType(0),
                                                                               mio::SimulationTime<double>(0.0),
                                                                               {0},
                                                                               Eigen::VectorXd::Constant(2, 1.0)},
                                                  mio::DampingSampling<double>{mio::UncertainValue<double>(0.25),
                                                                               mio::DampingLevel(1),
                                                                               mio::DampingType(0),
                                                                               mio::SimulationTime<double>(1.0),
                                                                               {
                                                                                   0,
                                                                                   1,
                                                                               },
                                                                               Eigen::VectorXd::Constant(2, 1.0)}};
    auto cmg = mio::ContactMatrixGroup<double>(2, 2);

    mio::apply_dampings(cmg, ds, [](auto&& v) {
        return mio::make_contact_damping_matrix(v);
    });

    ASSERT_THAT(
        cmg[0].get_dampings(),
        testing::ElementsAre(mio::SquareDamping<double>(Eigen::MatrixXd::Constant(2, 2, 0.5), mio::DampingLevel(0),
                                                        mio::DampingType(0), mio::SimulationTime<double>(0.0)),
                             mio::SquareDamping<double>(Eigen::MatrixXd::Constant(2, 2, 0.25), mio::DampingLevel(1),
                                                        mio::DampingType(0), mio::SimulationTime<double>(1.0))));
    ASSERT_THAT(cmg[1].get_dampings(), testing::ElementsAre(mio::SquareDamping<double>(
                                           Eigen::MatrixXd::Constant(2, 2, 0.25), mio::DampingLevel(1),
                                           mio::DampingType(0), mio::SimulationTime<double>(1.0))));
}

TEST(TestDampingSampling, contactMask)
{
    auto m = mio::make_contact_damping_matrix((Eigen::VectorXd(2) << 0.0, 0.5).finished()).eval();
    ASSERT_THAT(print_wrap(m),
                MatrixNear((Eigen::MatrixXd(2, 2) << 0.0, 1 - sqrt(0.5), 1 - sqrt(0.5), 0.5).finished()));
}

TEST(TestDampingSampling, mobilityMask)
{
    auto m = mio::make_mobility_damping_vector(mio::ColumnVectorShape<double>(6),
                                               (Eigen::VectorXd(2) << 0.5, 0.25).finished())
                 .eval();
    ASSERT_THAT(print_wrap(m), MatrixNear((Eigen::VectorXd(6) << 0.5, 0.5, 0.5, 0.25, 0.25, 0.25).finished()));
}
