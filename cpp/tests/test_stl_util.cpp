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
#include "memilio/utils/stl_util.h"
#include "memilio/utils/compiler_diagnostics.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <algorithm>

TEST(TestRange, index_operator)
{
    auto v = std::vector<int>{0, 1, 2, 3, 4, 5, 6};
    auto r = mio::Range(begin(v), end(v));

    EXPECT_EQ(v.size(), r.size());
    for (size_t i = 0; i < r.size(); ++i) {
        EXPECT_EQ(v[i], r[i]);
    }
}

TEST(TestRange, iterators)
{
    auto v = std::vector<int>{0, 1, 2, 3, 4, 5, 6};
    auto r = mio::Range(begin(v), end(v));

    EXPECT_THAT(r, testing::ElementsAreArray(v));
}

TEST(TestRange, reverse_iterators)
{
    auto v  = std::vector<int>{0, 1, 2, 3, 4, 5};
    auto r  = mio::Range(begin(v), end(v));
    auto v2 = std::vector<int>(r.rbegin(), r.rend());

    EXPECT_THAT(v2, testing::ElementsAre(5, 4, 3, 2, 1, 0));
}

TEST(TestRange, c_array)
{
    int v[] = {1, 2, 3, 4, 5, 6};
    auto r  = mio::Range(std::begin(v), std::end(v));

    EXPECT_THAT(r, testing::ElementsAreArray(v));
}

TEST(TestRange, reference_semantics)
{
    auto v = std::vector<int>{3, 4, 1, 2, 6, 7};
    auto r = mio::Range(begin(v), end(v));
    std::ranges::sort(v);

    EXPECT_THAT(r, testing::ElementsAreArray(v));
}

TEST(TestRange, partial_view)
{
    auto v = std::vector<int>{0, 1, 2, 3, 4, 5, 6};
    auto r = mio::Range(begin(v) + 2, end(v) - 1);

    EXPECT_THAT(r, testing::ElementsAre(2, 3, 4, 5));
}

namespace
{
struct Foo {
};

GCC_CLANG_DIAGNOSTIC(push)
GCC_CLANG_DIAGNOSTIC(ignored "-Wunneeded-internal-declaration")
std::ostream& operator<<(std::ostream& os, const Foo&)
{
    return os;
}
GCC_CLANG_DIAGNOSTIC(pop)

struct Bar {
};
} // namespace

TEST(TestTemplateUtils, has_stream_op)
{
    EXPECT_TRUE(mio::HasOstreamOperator<Foo>);
    EXPECT_FALSE(mio::HasOstreamOperator<Bar>);
}

TEST(TestInsertSortedReplace, normal)
{
    std::vector<int> v = {5};
    mio::insert_sorted_replace(v, 1);
    mio::insert_sorted_replace(v, 7);
    mio::insert_sorted_replace(v, 6);
    mio::insert_sorted_replace(v, 2);

    EXPECT_THAT(v, testing::ElementsAre(1, 2, 5, 6, 7));
}

TEST(TestInsertSortedReplace, reverse)
{
    std::vector<int> v = {5};
    auto pred          = [](auto&& l, auto r) {
        return r < l;
    };
    mio::insert_sorted_replace(v, 1, pred);
    mio::insert_sorted_replace(v, 7, pred);
    mio::insert_sorted_replace(v, 6, pred);
    mio::insert_sorted_replace(v, 2, pred);

    EXPECT_THAT(v, testing::ElementsAre(7, 6, 5, 2, 1));
}

TEST(TestInsertSortedReplace, replace)
{
    struct Foo_ {
        Foo_(int i1_, int i2_)
            : i1(i1_)
            , i2(i2_)
        {
        }
        int i1, i2;
        bool operator==(const Foo_& b) const
        {
            return i1 == b.i1 && i2 == b.i2;
        }
    };

    std::vector<Foo_> v = {{1, 1}, {3, 1}};
    auto pred           = [](auto&& l, auto r) {
        return l.i1 < r.i1;
    };
    mio::insert_sorted_replace(v, {2, 1}, pred);
    mio::insert_sorted_replace(v, {2, 2}, pred);
    mio::insert_sorted_replace(v, {1, 2}, pred);

    EXPECT_THAT(v, testing::ElementsAre(Foo_(1, 2), Foo_(2, 2), Foo_(3, 1)));
}

TEST(TestPathJoin, joinOne)
{
    EXPECT_EQ(mio::path_join("."), ".");
}

TEST(TestPathJoin, joinTwoMixedClasses)
{
    EXPECT_EQ(mio::path_join(".", "dir"), "./dir");
    EXPECT_EQ(mio::path_join("./", std::string("dir")), "./dir");
    EXPECT_EQ(mio::path_join(std::string("/"), "dir"), "/dir");
    EXPECT_EQ(mio::path_join(std::string("."), std::string("dir")), "./dir");
}

TEST(TestPathJoin, ignoreEmpty)
{
    EXPECT_EQ(mio::path_join(""), "");
    EXPECT_EQ(mio::path_join("", "dir"), "dir");
    EXPECT_EQ(mio::path_join("", "", "dir"), "dir");
    EXPECT_EQ(mio::path_join(".", "", "", "dir"), "./dir");
    EXPECT_EQ(mio::path_join("./", "", "", "dir"), "./dir");
}

namespace
{
class Base
{
public:
    virtual ~Base()
    {
    }
};

class Derived : public Base
{
};
} // namespace

TEST(TestDynamicUniquePtrCast, notNull)
{
    std::unique_ptr<Base> upb = std::make_unique<Derived>();
    auto pb                   = upb.get();
    auto upd                  = mio::dynamic_unique_ptr_cast<Derived>(std::move(upb));
    EXPECT_EQ(upd.get(), dynamic_cast<Derived*>(pb));
}

TEST(TestDynamicUniquePtrCast, null)
{
    std::unique_ptr<Base> upb;
    auto upd = mio::dynamic_unique_ptr_cast<Derived>(std::move(upb));
    EXPECT_EQ(upd.get(), nullptr);
}

TEST(TestDynamicUniquePtrCast, notDerived)
{
    std::unique_ptr<Base> upb = std::make_unique<Base>();
    auto upd                  = mio::dynamic_unique_ptr_cast<Derived>(std::move(upb));
    EXPECT_EQ(upd.get(), nullptr);
}

TEST(TestContains, normalCase)
{
    auto v = std::vector<int>{4, 1, 3, 6, 9};
    ASSERT_TRUE(mio::contains(v.begin(), v.end(), [](auto&& e) {
        return e == 3;
    }));
    ASSERT_FALSE(mio::contains(v.begin(), v.end(), [](auto&& e) {
        return e == 7;
    }));
}

TEST(TestContains, empty)
{
    auto v = std::vector<int>();
    ASSERT_FALSE(mio::contains(v.begin(), v.end(), [](auto&&) {
        return true;
    }));
}

TEST(EnumMembers, works)
{
    enum class E
    {
        A,
        B,
        Count
    };
    ASSERT_THAT(mio::enum_members<E>(), testing::ElementsAre(E::A, E::B));
}

TEST(TestContains, set_ostream_format)
{
    std::ostringstream output;

    mio::set_ostream_format(output, 10, 2, '*');
    output << 3.14159;
    std::string expected_output_1 = "******3.14";
    std::string actual_output_1   = output.str();
    EXPECT_EQ(expected_output_1, actual_output_1);

    output.str("");

    mio::set_ostream_format(output, 7, 3, '#');
    output << 42.12345678;
    std::string expected_output_2 = "#42.123";
    std::string actual_output_2   = output.str();
    EXPECT_EQ(expected_output_2, actual_output_2);

    output.str("");

    mio::set_ostream_format(output, 8, 4);
    output << 123.456;
    std::string expected_output_3 = "123.4560";
    std::string actual_output_3   = output.str();
    EXPECT_EQ(expected_output_3, actual_output_3);
}
