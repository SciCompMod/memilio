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
#include "memilio/utils/random_number_generator.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <cstdint>
#include <iterator>

TEST(RandomNumberGenerator, counter)
{
    mio::RandomNumberGenerator rng;
    uint64_t n = 10;
    for (auto i = uint64_t(0); i < n; ++i) {
        rng();
    }
    ASSERT_EQ(rng.get_counter(), mio::Counter<uint64_t>(n)); //counter was incremented

    auto s1 = rng();
    for (auto i = uint64_t(0); i < n; ++i) {
        rng();
    }
    ASSERT_EQ(rng.get_counter(), mio::Counter<uint64_t>(2 * n + 1)); //counter was further incremented

    rng.set_counter(mio::Counter<uint64_t>{n});
    ASSERT_EQ(rng.get_counter(), mio::Counter<uint64_t>(n)); //counter was reverted
    auto s2 = rng();

    ASSERT_EQ(s1, s2);
}

TEST(RandomNumberGenerator, subsequences)
{
    uint32_t subsequence_index{5};
    mio::Counter<uint32_t> subsequence_counter{12};

    auto counter = mio::rng_totalsequence_counter<uint64_t>(subsequence_index, subsequence_counter);
    ASSERT_EQ(counter, mio::Counter<uint64_t>{5 * 4'294'967'296ll + 12});

    auto subsequence_counter2 = mio::rng_subsequence_counter<uint32_t>(counter);
    ASSERT_EQ(subsequence_counter2, subsequence_counter);
}

TEST(TestDiscreteDistribution, generate)
{
    auto distribution = mio::DiscreteDistributionInPlace<size_t>();
    auto rng          = mio::RandomNumberGenerator();

    std::vector<double> weights;
    for (size_t i = 0; i < 50; i++) {
        weights = {};
        ASSERT_EQ(distribution(rng, {weights}), 0);

        weights = {0.5};
        ASSERT_EQ(distribution(rng, {weights}), 0);

        weights = {0.5, 1.3, 0.1, 0.4, 0.3};
        auto d  = distribution(rng, {weights});
        ASSERT_GE(d, 0);
        ASSERT_LE(d, 4);
    }
}
