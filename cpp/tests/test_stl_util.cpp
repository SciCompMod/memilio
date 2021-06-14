#include <epidemiology/utils/stl_util.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

TEST(TestRange, index_operator)
{
    auto v = std::vector<int>{0, 1, 2, 3, 4, 5, 6};
    auto r = epi::make_range(begin(v), end(v));

    EXPECT_EQ(v.size(), r.size());
    for (size_t i = 0; i < r.size(); ++i) {
        EXPECT_EQ(v[i], r[i]);
    }
}

TEST(TestRange, iterators)
{
    auto v = std::vector<int>{0, 1, 2, 3, 4, 5, 6};
    auto r = epi::make_range(begin(v), end(v));

    EXPECT_THAT(r, testing::ElementsAreArray(v));
}

TEST(TestRange, reverse_iterators)
{
    auto v = std::vector<int>{0, 1, 2, 3, 4, 5};
    auto r = epi::make_range(begin(v), end(v));
    auto v2 = std::vector<int>(r.rbegin(), r.rend());

    EXPECT_THAT(v2, testing::ElementsAre(5, 4, 3, 2, 1, 0));
}

TEST(TestRange, c_array)
{
    int v[] = {1, 2, 3, 4, 5, 6};
    auto r  = epi::make_range(std::begin(v), std::end(v));

    EXPECT_THAT(r, testing::ElementsAreArray(v));
}

TEST(TestRange, reference_semantics)
{
    auto v = std::vector<int>{3, 4, 1, 2, 6, 7};
    auto r = epi::make_range(begin(v), end(v));
    std::sort(begin(v), end(v));

    EXPECT_THAT(r, testing::ElementsAreArray(v));
}

TEST(TestRange, partial_view)
{
    auto v = std::vector<int>{0, 1, 2, 3, 4, 5, 6};
    auto r = epi::make_range(begin(v) + 2, end(v) - 1);

    EXPECT_THAT(r, testing::ElementsAre(2, 3, 4, 5));
}

namespace
{
struct Foo {
};

std::ostream& operator<<(std::ostream& os, const Foo&)
{
    return os;
}

struct Bar {
};
} // namespace

TEST(TestTemplateUtils, has_stream_op)
{
    EXPECT_TRUE(epi::has_ostream_op<Foo>::value);
    EXPECT_FALSE(epi::has_ostream_op<Bar>::value);
}

TEST(TestInsertSortedReplace, normal)
{
    std::vector<int> v = {5};
    epi::insert_sorted_replace(v, 1);
    epi::insert_sorted_replace(v, 7);
    epi::insert_sorted_replace(v, 6);
    epi::insert_sorted_replace(v, 2);

    EXPECT_THAT(v, testing::ElementsAre(1, 2, 5, 6, 7));
}

TEST(TestInsertSortedReplace, returnsValidIterator)
{
    std::vector<int> v;
    int x;

    //There is no GTEST_NO_DEATH macro so we just let the test crash.
    //If this test crashes, the function does not return a valid iterator.
    //Dereferencing an invalid iterator is undefined behavior so the test 
    //may behave unexpectedly (pass, fail, or something else) if the iterator is invalid.
    x = *epi::insert_sorted_replace(v, 5);
    x = *epi::insert_sorted_replace(v, 1);
    x = *epi::insert_sorted_replace(v, 4);
    x = *epi::insert_sorted_replace(v, 7);
    ASSERT_EQ(x, 7);
}

TEST(TestInsertSortedReplace, reverse)
{
    std::vector<int> v = {5};
    auto pred          = [](auto&& l, auto r) {
        return r < l;
    };
    epi::insert_sorted_replace(v, 1, pred);
    epi::insert_sorted_replace(v, 7, pred);
    epi::insert_sorted_replace(v, 6, pred);
    epi::insert_sorted_replace(v, 2, pred);

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
    epi::insert_sorted_replace(v, {2, 1}, pred);
    epi::insert_sorted_replace(v, {2, 2}, pred);
    epi::insert_sorted_replace(v, {1, 2}, pred);

    EXPECT_THAT(v, testing::ElementsAre(Foo_(1, 2), Foo_(2, 2), Foo_(3, 1)));
}

TEST(TestPathJoin, joinOne)
{
    EXPECT_EQ(epi::path_join("."), ".");
}

TEST(TestPathJoin, joinTwoMixedClasses)
{
    EXPECT_EQ(epi::path_join(".", "dir"), "./dir");
    EXPECT_EQ(epi::path_join("./", std::string("dir")), "./dir");
    EXPECT_EQ(epi::path_join(std::string("/"), "dir"), "/dir");
    EXPECT_EQ(epi::path_join(std::string("."), std::string("dir")), "./dir");
}

TEST(TestPathJoin, ignoreEmpty)
{
    EXPECT_EQ(epi::path_join(""), "");
    EXPECT_EQ(epi::path_join("", "dir"), "dir");
    EXPECT_EQ(epi::path_join("", "", "dir"), "dir");
    EXPECT_EQ(epi::path_join(".", "", "", "dir"), "./dir");
    EXPECT_EQ(epi::path_join("./", "", "", "dir"), "./dir");
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
    auto upd                  = epi::dynamic_unique_ptr_cast<Derived>(std::move(upb));
    EXPECT_EQ(upd.get(), dynamic_cast<Derived*>(pb));
}

TEST(TestDynamicUniquePtrCast, null)
{
    std::unique_ptr<Base> upb;
    auto upd = epi::dynamic_unique_ptr_cast<Derived>(std::move(upb));
    EXPECT_EQ(upd.get(), nullptr);
}

TEST(TestDynamicUniquePtrCast, notDerived)
{
    std::unique_ptr<Base> upb = std::make_unique<Base>();
    auto upd                  = epi::dynamic_unique_ptr_cast<Derived>(std::move(upb));
    EXPECT_EQ(upd.get(), nullptr);
}

TEST(TestContains, normalCase)
{
    auto v = std::vector<int>{4, 1, 3, 6, 9};
    ASSERT_TRUE(epi::contains(v.begin(), v.end(), [](auto&& e) {
        return e == 3;
    }));
    ASSERT_FALSE(epi::contains(v.begin(), v.end(), [](auto&& e) {
        return e == 7;
    }));
}

TEST(TestContains, empty)
{
    auto v = std::vector<int>();
    ASSERT_FALSE(epi::contains(v.begin(), v.end(), [](auto&&) {
        return true;
    }));
}