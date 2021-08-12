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
#include "epidemiology/secir/damping_sampling.h"
#include "epidemiology/secir/contact_matrix.h"
#include "matchers.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <cmath>

TEST(TestDampingSampling, apply)
{
    auto ds  = std::vector<epi::DampingSampling>{epi::DampingSampling{epi::UncertainValue(0.5),
                                                                     epi::DampingLevel(0),
                                                                     epi::DampingType(0),
                                                                     epi::SimulationTime(0.0),
                                                                     {0},
                                                                     Eigen::VectorXd::Constant(2, 1.0)},
                                                epi::DampingSampling{epi::UncertainValue(0.25),
                                                                     epi::DampingLevel(1),
                                                                     epi::DampingType(0),
                                                                     epi::SimulationTime(1.0),
                                                                     {
                                                                         0,
                                                                         1,
                                                                     },
                                                                     Eigen::VectorXd::Constant(2, 1.0)}};
    auto cmg = epi::ContactMatrixGroup(2, 2);

    epi::apply_dampings(cmg, ds, [](auto&& v) {
        return epi::make_contact_damping_matrix(v);
    });

    ASSERT_THAT(cmg[0].get_dampings(),
                testing::ElementsAre(epi::SquareDamping(Eigen::MatrixXd::Constant(2, 2, 0.5), epi::DampingLevel(0),
                                                        epi::DampingType(0), epi::SimulationTime(0.0)),
                                     epi::SquareDamping(Eigen::MatrixXd::Constant(2, 2, 0.25), epi::DampingLevel(1),
                                                        epi::DampingType(0), epi::SimulationTime(1.0))));
    ASSERT_THAT(cmg[1].get_dampings(),
                testing::ElementsAre(epi::SquareDamping(Eigen::MatrixXd::Constant(2, 2, 0.25), epi::DampingLevel(1),
                                                        epi::DampingType(0), epi::SimulationTime(1.0))));
}

TEST(TestDampingSampling, contactMask)
{
    auto m = epi::make_contact_damping_matrix((Eigen::VectorXd(2) << 0.0, 0.5).finished()).eval();
    ASSERT_THAT(print_wrap(m), MatrixNear((Eigen::MatrixXd(2, 2) << 0.0, 1-sqrt(0.5), 1-sqrt(0.5), 0.5).finished()));
}

TEST(TestDampingSampling, migrationMask)
{
    auto m = epi::make_migration_damping_vector(epi::ColumnVectorShape(6),
                                                       (Eigen::VectorXd(2) << 0.5, 0.25).finished())
                 .eval();
    ASSERT_THAT(print_wrap(m), MatrixNear((Eigen::VectorXd(6) << 0.5, 0.5, 0.5, 0.25, 0.25, 0.25).finished()));
}