#include "memilio/utils/random_number_generator.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <algorithm>
#include <iterator>

TEST(RandomNumberGenerator, blocks)
{
    mio::RandomNumberGenerator rng1, rng2;
    std::vector<unsigned int> seeds = { 1, 2, 3, 4 };
    rng1.seed(seeds);
    rng2.seed(seeds);

    std::vector<mio::RandomNumberGenerator::result_type> samples1;
    std::generate_n(std::back_inserter(samples1), 100, [&rng1] { return rng1(); });

    rng2.set_block_size(5);
    std::vector<mio::RandomNumberGenerator::result_type> dump;
    std::generate_n(std::back_inserter(dump), 23, [&rng2] { return rng2(); });
    rng2.forward_to_block(19);
    std::vector<mio::RandomNumberGenerator::result_type> samples2;
    std::generate_n(std::back_inserter(samples2), 5, [&rng2] { return rng2(); });

    ASSERT_THAT(samples2, testing::ElementsAreArray(samples1.begin() + 95, samples1.end()));
}

#ifndef NDEBUG
TEST(RandomNumberGenerator, block_error)
{
    mio::RandomNumberGenerator rng;
    rng.set_block_size(5);
    std::vector<mio::RandomNumberGenerator::result_type> dump;
    std::generate_n(std::back_inserter(dump), 123, [&rng] {
        return rng();
    });
    ASSERT_DEATH(rng.forward_to_block(19), ".*");
}
#endif
