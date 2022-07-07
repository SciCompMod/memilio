#include "memilio/utils/transform_iterator.h"
#include "memilio/utils/stl_util.h"
#include "gtest/gtest.h"
#include "gmock/gmock.h"

TEST(TestTransformIterator, produce_lvalue)
{
    struct Foo
    {
        int& get_i() {
             return i;
        }
        int i;
    };
    auto v = std::vector<Foo>{{1}, {2}};
    auto get_i = [](auto&& foo) -> auto& {
        return foo.get_i();
    };
    auto b = mio::make_transform_iterator(v.begin(), get_i);
    auto e = mio::make_transform_iterator(v.end(), get_i);
    (*(b + 1)) = 3; //element can be modified through the iterator
    ASSERT_THAT(mio::make_range(b, e), testing::ElementsAre(1, 3));
}

TEST(TestTransformIterator, produce_rvalue)
{
    auto v = std::vector<int>{1, 2};
    auto square = [](auto&& i) {
        return i * i;
    };
    auto vt = std::vector<int>(mio::make_transform_iterator(v.begin(), square), mio::make_transform_iterator(v.end(), square));
    ASSERT_THAT(vt, testing::ElementsAre(1, 4));
}
