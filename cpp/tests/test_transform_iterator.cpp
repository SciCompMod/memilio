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
#include "memilio/utils/transform_iterator.h"
#include "memilio/utils/stl_util.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

TEST(TestTransformIterator, produce_lvalue)
{
    struct Foo {
        int& get_i()
        {
            return i;
        }
        int i;
    };
    auto v     = std::vector<Foo>{{1}, {2}};
    auto get_i = [](auto&& foo) -> auto& {
        return foo.get_i();
    };
    auto b     = mio::make_transform_iterator(v.begin(), get_i);
    auto e     = mio::make_transform_iterator(v.end(), get_i);
    (*(b + 1)) = 3; //element can be modified through the iterator
    ASSERT_THAT(mio::make_range(b, e), testing::ElementsAre(1, 3));
}

TEST(TestTransformIterator, produce_rvalue)
{
    auto v      = std::vector<int>{1, 2};
    auto square = [](auto&& i) {
        return i * i;
    };
    auto vt = std::vector<int>(mio::make_transform_iterator(v.begin(), square),
                               mio::make_transform_iterator(v.end(), square));
    ASSERT_THAT(vt, testing::ElementsAre(1, 4));
}
