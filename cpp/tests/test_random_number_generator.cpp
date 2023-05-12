#include "memilio/utils/random_number_generator.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <iterator>

TEST(RandomNumberGenerator, set_counter)
{
    mio::RandomNumberGenerator rng;
    uint64_t n = 10;
    for (auto i = uint64_t(0); i < n; ++i) {
        rng();
    }
    ASSERT_EQ(rng.get_counter(), mio::Counter<uint64_t>(n));
    auto s1 = rng();
    for (auto i = uint64_t(0); i < n; ++i) {
        rng();
    }
    ASSERT_EQ(rng.get_counter(), mio::Counter<uint64_t>(2 * n + 1));

    rng.set_counter(mio::Counter<uint64_t>{n});
    ASSERT_EQ(rng.get_counter(), mio::Counter<uint64_t>(n));
    auto s2 = rng();

    ASSERT_EQ(s1, s2);
}

TEST(TestDiscreteDistribution, generate)
{
    auto distribution = mio::DiscreteDistributionInPlace<size_t>();
    auto rng = mio::RandomNumberGenerator();

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
