/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Julia Bicker
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
#include "abm/personal_rng.h"
#include "memilio/utils/abstract_parameter_distribution.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/random_number_generator.h"
#include "models/abm/personal_rng.h"
#include <gtest/gtest.h>
#include <vector>

TEST(AbstractParameterDist, test_abstract_normal_distribution)
{
    mio::AbstractParameterDistribution p1(mio::ParameterDistributionNormal(10., 0.5));
    mio::AbstractParameterDistribution p2(mio::ParameterDistributionNormal(20., 0.3));
    auto params1 = std::vector<double>{10., 0.5};
    auto params2 = std::vector<double>{20., 0.3};
    //Check parameters
    EXPECT_EQ(p1.params(), params1);
    EXPECT_EQ(p2.params(), params2);
    //check sampling function
    auto counter = mio::Counter<uint32_t>{0};
    auto prng    = mio::abm::PersonalRandomNumberGenerator(mio::thread_local_rng().get_key(), 0, counter);
    EXPECT_GE(p2.get(mio::thread_local_rng()), 0.);
    EXPECT_GE(p2.get(prng), 0.);
    //check smaller impl
    mio::AbstractParameterDistribution p3(mio::ParameterDistributionUniform(0., 1.));
    bool smaller     = (p1 < p2);
    bool not_smaller = (p1 < p3);
    EXPECT_EQ(smaller, true);
    EXPECT_EQ(not_smaller, false);
}

TEST(AbstractParameterDist, test_abstract_uniform_distribution)
{
    mio::AbstractParameterDistribution p1(mio::ParameterDistributionUniform(0., 1.));
    mio::AbstractParameterDistribution p2(mio::ParameterDistributionUniform(2., 3.));
    auto params1 = std::vector<double>{0., 1.};
    auto params2 = std::vector<double>{2., 3.};
    //Check parameters
    EXPECT_EQ(p1.params(), params1);
    EXPECT_EQ(p2.params(), params2);
    //check sampling function
    auto counter = mio::Counter<uint32_t>{0};
    auto prng    = mio::abm::PersonalRandomNumberGenerator(mio::thread_local_rng().get_key(), 0, counter);
    EXPECT_GE(p1.get(mio::thread_local_rng()), 0.);
    EXPECT_GE(p1.get(prng), 0.);
    EXPECT_LE(p1.get(mio::thread_local_rng()), 1.);
    EXPECT_LE(p1.get(prng), 1.);
    //check smaller impl
    mio::AbstractParameterDistribution p3(mio::ParameterDistributionNormal(0., 1.));
    bool smaller     = (p1 < p2);
    bool not_smaller = (p1 < p3);
    EXPECT_EQ(smaller, true);
    EXPECT_EQ(not_smaller, false);
}

TEST(AbstractParameterDist, test_abstract_lognormal_distribution)
{
    mio::AbstractParameterDistribution p1(mio::ParameterDistributionLogNormal(1., 0.25));
    mio::AbstractParameterDistribution p2(mio::ParameterDistributionLogNormal(2., 0.1));
    auto params1 = std::vector<double>{1., 0.25};
    auto params2 = std::vector<double>{2., 0.1};
    //Check parameters
    EXPECT_EQ(p1.params(), params1);
    EXPECT_EQ(p2.params(), params2);
    //check sampling function
    auto counter = mio::Counter<uint32_t>{0};
    auto prng    = mio::abm::PersonalRandomNumberGenerator(mio::thread_local_rng().get_key(), 0, counter);
    EXPECT_GE(p1.get(mio::thread_local_rng()), 0.);
    EXPECT_GE(p1.get(prng), 0.);
    //check smaller impl
    mio::AbstractParameterDistribution p3(mio::ParameterDistributionNormal(0., 1.));
    bool smaller     = (p1 < p2);
    bool not_smaller = (p1 < p3);
    EXPECT_EQ(smaller, true);
    EXPECT_EQ(not_smaller, false);
}

TEST(AbstractParameterDist, test_abstract_exponential_distribution)
{
    mio::AbstractParameterDistribution p1(mio::ParameterDistributionExponential(1.));
    mio::AbstractParameterDistribution p2(mio::ParameterDistributionExponential(2.));
    auto params1 = std::vector<double>{1.};
    auto params2 = std::vector<double>{2.};
    //Check parameters
    EXPECT_EQ(p1.params(), params1);
    EXPECT_EQ(p2.params(), params2);
    //check sampling function
    auto counter = mio::Counter<uint32_t>{0};
    auto prng    = mio::abm::PersonalRandomNumberGenerator(mio::thread_local_rng().get_key(), 0, counter);
    EXPECT_GE(p1.get(mio::thread_local_rng()), 0.);
    EXPECT_GE(p1.get(prng), 0.);
    //check smaller impl
    mio::AbstractParameterDistribution p3(mio::ParameterDistributionNormal(0., 1.));
    bool smaller     = (p1 < p2);
    bool not_smaller = (p1 < p3);
    EXPECT_EQ(smaller, true);
    EXPECT_EQ(not_smaller, false);
}

TEST(AbstractParameterDist, test_abstract_constant_distribution)
{
    mio::AbstractParameterDistribution p1(mio::ParameterDistributionConstant(1.));
    mio::AbstractParameterDistribution p2(mio::ParameterDistributionConstant(2.));
    auto params1 = std::vector<double>{1.};
    auto params2 = std::vector<double>{2.};
    //Check parameters
    EXPECT_EQ(p1.params(), params1);
    EXPECT_EQ(p2.params(), params2);
    //check sampling function
    auto counter = mio::Counter<uint32_t>{0};
    auto prng    = mio::abm::PersonalRandomNumberGenerator(mio::thread_local_rng().get_key(), 0, counter);
    EXPECT_EQ(p1.get(mio::thread_local_rng()), 1.);
    EXPECT_EQ(p1.get(prng), 1.);
    //check smaller impl
    mio::AbstractParameterDistribution p3(mio::ParameterDistributionNormal(0., 1.));
    bool smaller     = (p1 < p2);
    bool not_smaller = (p1 < p3);
    EXPECT_EQ(smaller, true);
    EXPECT_EQ(not_smaller, false);
}
