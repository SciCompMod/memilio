/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: René Schmieding
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
#include "load_test_data.h"
#include "memilio/io/json_serializer.h"
#include "memilio/utils/random_number_generator.h"
#include <gtest/gtest.h>

template <class TestType>
class TestRNG : public ::testing::Test
{
protected:
    TestRNG();

public:
    const std::vector<unsigned> seeds{1214370497, 2300667602, 3141592653, 2158220490, 1234567890, 3429711942};
    std::vector<typename TestType::result_type> compare; // drawn samples in random_number_generator.csv
    TestType rng;

private:
    std::vector<typename TestType::result_type> load_random_numbers()
    {
        std::vector<typename TestType::result_type> comp;
        std::ifstream result_file(get_test_data_file_path("random_numbers.csv"));
        std::string reader;
        getline(result_file, reader);
        std::stringstream parser(reader);
        typename TestType::result_type r;
        while (parser >> r) {
            comp.push_back(r);
        }
        result_file.close();
        return comp;
    }
};

template <>
TestRNG<mio::RandomNumberGenerator>::TestRNG()
    : compare(load_random_numbers())
{
    rng.seed(seeds);
}

template <>
TestRNG<mio::CachedRNG>::TestRNG()
    : compare(load_random_numbers())
    , rng(compare.size(), std::vector<mio::CachedRNG::result_type>{}, seeds)
{
}

using TestTypes = ::testing::Types<mio::RandomNumberGenerator, mio::CachedRNG>;

TYPED_TEST_SUITE(TestRNG, TestTypes);

TYPED_TEST(TestRNG, Bounds)
{
    EXPECT_LE(this->rng.min(), this->rng.max());
}

TYPED_TEST(TestRNG, Seeding)
{
    using result_type = typename TypeParam::result_type;
    std::vector<result_type> samples(this->compare.size());
    // compare seeds
    auto s = this->rng.get_seeds();
    EXPECT_EQ(s.size(), this->seeds.size());
    for (unsigned i = 0; i < s.size(); i++) {
        EXPECT_EQ(s[i], this->seeds[i]) << "i = " << i;
    }
    // draw and compare random numbers
    for (auto& x : samples) {
        x = this->rng();
    }
    for (unsigned i = 0; i < this->compare.size(); i++) {
        EXPECT_EQ(samples[i], this->compare[i]) << "i = " << i;
    }
}
