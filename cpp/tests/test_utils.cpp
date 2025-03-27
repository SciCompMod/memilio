/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding 
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
#include "memilio/utils/index.h"
#include "memilio/utils/index_range.h"

#include <gtest/gtest.h>

template <size_t Tag>
struct CategoryTag : public mio::Index<CategoryTag<Tag>> {
    CategoryTag(size_t value)
        : mio::Index<CategoryTag<Tag>>(value)
    {
    }
};

TEST(TestUtils, reduce_index)
{
    // the superindex
    mio::Index<CategoryTag<1>, CategoryTag<2>, CategoryTag<3>> i{CategoryTag<1>(1), CategoryTag<2>(2),
                                                                 CategoryTag<3>(3)};
    // some subindices
    mio::Index<CategoryTag<1>, CategoryTag<2>> reference_i12{CategoryTag<1>(1), CategoryTag<2>(2)};
    mio::Index<CategoryTag<2>, CategoryTag<3>> reference_i23{CategoryTag<2>(2), CategoryTag<3>(3)};
    mio::Index<CategoryTag<1>, CategoryTag<3>> reference_i13{CategoryTag<1>(1), CategoryTag<3>(3)};
    // test reduction results
    auto result_i12 = mio::reduce_index<decltype(reference_i12)>(i);
    EXPECT_EQ(result_i12, reference_i12);
    auto result_i23 = mio::reduce_index<decltype(reference_i23)>(i);
    EXPECT_EQ(result_i23, reference_i23);
    auto result_i13 = mio::reduce_index<decltype(reference_i13)>(i);
    EXPECT_EQ(result_i13, reference_i13);

    // test reduction of an index with non-unique categories
    mio::Index<CategoryTag<1>, CategoryTag<1>, CategoryTag<1>> j{CategoryTag<1>(1), CategoryTag<1>(2),
                                                                 CategoryTag<1>(3)};
    mio::Index<CategoryTag<1>, CategoryTag<1>> reference_j{CategoryTag<1>(1), CategoryTag<1>(1)};
    auto result_j = mio::reduce_index<decltype(reference_j)>(j);
    EXPECT_EQ(result_j, reference_j);
}

TEST(TestUtils, extend_index)
{
    // the superindex
    mio::Index<CategoryTag<1>, CategoryTag<2>, CategoryTag<3>> reference_i{CategoryTag<1>(1), CategoryTag<2>(2),
                                                                           CategoryTag<3>(3)};
    // some subindices
    mio::Index<CategoryTag<1>, CategoryTag<2>> i12{CategoryTag<1>(1), CategoryTag<2>(2)};
    mio::Index<CategoryTag<2>, CategoryTag<3>> i23{CategoryTag<2>(2), CategoryTag<3>(3)};
    mio::Index<CategoryTag<1>, CategoryTag<3>> i13{CategoryTag<1>(1), CategoryTag<3>(3)};
    // test extension results
    auto result_i12 = mio::extend_index<decltype(reference_i)>(i12, 3);
    EXPECT_EQ(result_i12, reference_i);
    auto result_i23 = mio::extend_index<decltype(reference_i)>(i23, 1);
    EXPECT_EQ(result_i23, reference_i);
    auto result_i13 = mio::extend_index<decltype(reference_i)>(i13, 2);
    EXPECT_EQ(result_i13, reference_i);

    // test extension of an index with non-unique categories
    mio::Index<CategoryTag<1>, CategoryTag<1>, CategoryTag<1>> j{CategoryTag<1>(1), CategoryTag<1>(2),
                                                                 CategoryTag<1>(3)};
    mio::Index<CategoryTag<1>, CategoryTag<1>, CategoryTag<1>, CategoryTag<1>> reference_j{
        CategoryTag<1>(1), CategoryTag<1>(1), CategoryTag<1>(1), CategoryTag<1>(1)};
    auto result_j = mio::extend_index<decltype(reference_j)>(j);
    EXPECT_EQ(result_j, reference_j);
}

TEST(TestUtils, IndexRange)
{
    using I = mio::Index<CategoryTag<1>, CategoryTag<2>, CategoryTag<3>>;
    I dims{CategoryTag<1>(2), CategoryTag<2>(3), CategoryTag<3>(5)};
    mio::IndexRange<I> range(dims);

    I reference_begin{CategoryTag<1>(0), CategoryTag<2>(0), CategoryTag<3>(0)};
    I reference_end{CategoryTag<1>(2), CategoryTag<2>(0), CategoryTag<3>(0)};
    // check begin / end iterators
    EXPECT_EQ(*range.begin(), reference_begin);
    EXPECT_EQ(*range.end(), reference_end);
    // test increments
    auto iterator = range.begin();
    for (size_t reference_flat_index = 0; reference_flat_index < 2 * 3 * 5; reference_flat_index++) {
        // manually compute flatten_index(*iterator, dims)
        size_t flat_index =
            (3 * 5) * mio::get<0>(*iterator).get() + (5) * mio::get<1>(*iterator).get() + mio::get<2>(*iterator).get();
        EXPECT_EQ(flat_index, reference_flat_index);
        iterator++;
    }
}