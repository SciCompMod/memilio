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
#include "memilio/utils/memory.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/parameter_distributions.h"
#include <distributions_helpers.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <memory>

TEST(TestUncertain, uncertain_value_basic)
{
    mio::UncertainValue<double> val(3.0);
    EXPECT_EQ(val, 3.0);

    val = 2.0;
    EXPECT_EQ(val, 2.0);
}

TEST(TestUncertain, uncertain_value_copy)
{
    mio::UncertainValue<double> val(2.0);
    double dev_rel     = 0.2;
    double lower_bound = std::max(1e-6, (1 - dev_rel * 2.6) * val);
    double upper_bound = (1 + dev_rel * 2.6) * val;
    val.set_distribution(mio::ParameterDistributionNormal(lower_bound, upper_bound, val, dev_rel * val));

    mio::UncertainValue<double> val2(val);
    EXPECT_EQ(val2, 2.0);

    EXPECT_NE(val.get_distribution().get(), val2.get_distribution().get()); // dists get copied
    check_distribution(*val.get_distribution().get(), *val2.get_distribution().get());
}

TEST(TestUncertain, random_sample)
{
    mio::UncertainValue<double> val(2.0);

    auto mock_dist_ref = MockParameterDistributionRef<testing::StrictMock<MockParameterDistribution>>();
    EXPECT_CALL(mock_dist_ref.get_mock(), get_rand_sample())
        .Times(2)
        .WillOnce(testing::Return(0.5))
        .WillOnce(testing::Return(1.5));
    val.set_distribution(mock_dist_ref);

    EXPECT_EQ(val, 2.0);
    val.draw_sample();
    EXPECT_EQ(val, 0.5);
    val = 1.0;
    EXPECT_EQ(val, 1.0);
    val.draw_sample();
    EXPECT_EQ(val, 1.5);
}

TEST(TestUncertain, uncertain_value_assign)
{
    mio::UncertainValue<double> val(3.0);
    double dev_rel     = 0.2;
    double lower_bound = std::max(1e-6, (1 - dev_rel * 2.6) * val);
    double upper_bound = (1 + dev_rel * 2.6) * val;
    val.set_distribution(mio::ParameterDistributionNormal(lower_bound, upper_bound, val, dev_rel * val));

    double& dval = val;
    dval         = 4.0;
    EXPECT_EQ(val, 4.0);

    mio::UncertainValue<double> val2 = 4.0;
    EXPECT_EQ(val2, val); // only checks doubles, not dists

    mio::UncertainValue<double> val3;
    val3 = val;

    EXPECT_EQ(val3, val);
    check_distribution(*val3.get_distribution().get(), *val.get_distribution().get());
}

TEST(TestUncertain, uncertain_value_predef)
{
    mio::UncertainValue<double> val(3.0);
    double dev_rel     = 0.2;
    double lower_bound = std::max(1e-6, (1 - dev_rel * 2.6) * val);
    double upper_bound = (1 + dev_rel * 2.6) * val;
    val.set_distribution(mio::ParameterDistributionNormal(lower_bound, upper_bound, val, dev_rel * val));

    val.get_distribution().get()->add_predefined_sample(17111.0);
    mio::UncertainValue<double> val2(val);
    val.draw_sample();
    EXPECT_EQ(val, 17111.0);

    val.draw_sample();
    EXPECT_NE(val, 17111.0);

    val2.draw_sample();
    EXPECT_EQ(val2, 17111.0);
}

TEST(TestUncertain, uncertain_matrix)
{
    mio::ContactMatrix<double> contact_matrix(Eigen::MatrixXd::NullaryExpr(2, 2, [](auto i, auto j) -> double {
        return (i + 1) * (j + 1);
    }));
    contact_matrix.add_damping(0.7, mio::SimulationTime<double>(30.));

    mio::UncertainContactMatrix<double> uncertain_mat{{contact_matrix}};

    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_baseline()(0, 1), 2);
    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_baseline()(1, 1), 4);
    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_dampings()[0].get_coeffs()(1, 1), 0.7);

    uncertain_mat.get_dampings().emplace_back(mio::UncertainValue<double>(0.5), mio::DampingLevel(0),
                                              mio::DampingType(0), mio::SimulationTime<double>(3.0),
                                              std::vector<size_t>(1, size_t(0)), Eigen::VectorXd::Constant(2, 1.0));

    uncertain_mat.get_school_holiday_damping() = mio::DampingSampling<double>(
        mio::UncertainValue<double>(1.), mio::DampingLevel(1), mio::DampingType(0), mio::SimulationTime<double>(0.0),
        std::vector<size_t>(1, size_t(0)), (Eigen::VectorXd(2) << 1.0, 0.0).finished());
    uncertain_mat.get_school_holidays().assign({{mio::SimulationTime<double>(5.0), mio::SimulationTime<double>(17.0)}});

    uncertain_mat.draw_sample(true);
    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_dampings().size(), 4);

    uncertain_mat.draw_sample();
    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_dampings().size(), 3);

    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_dampings()[0].get_level(), mio::DampingLevel(0));
    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_dampings()[0].get_coeffs()(0, 0), 0.5);
    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_dampings()[0].get_time(), mio::SimulationTime<double>(3.0));

    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_dampings()[1].get_level(), mio::DampingLevel(1));
    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_dampings()[1].get_coeffs()(0, 0), 1.0);
    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_dampings()[1].get_time(), mio::SimulationTime<double>(5.0));

    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_dampings()[2].get_level(), mio::DampingLevel(1));
    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_dampings()[2].get_coeffs()(0, 0), 0.0);
    EXPECT_EQ(uncertain_mat.get_cont_freq_mat()[0].get_dampings()[2].get_time(), mio::SimulationTime<double>(17.0));
}
