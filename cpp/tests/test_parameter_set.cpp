#include <epidemiology/utils/parameter_set.h>
#include <gtest/gtest.h>
#include <gmock/gmock.h>

struct IntParam1 {
    using Type = int;
    static constexpr Type get_default()
    {
        return 1;
    }
};

struct IntParam2 {
    using Type = int;
    static constexpr Type get_default()
    {
        return 2;
    }
};

struct DoubleParam {
    using Type = double;
    static constexpr Type get_default()
    {
        return 1.0;
    }
};

struct StringParam {
    using Type = std::string;
};

struct NotDefaultConstructible {
    int i;
    NotDefaultConstructible(int i_)
        : i(i_)
    {
    }
    ~NotDefaultConstructible()
    {
    }
    bool operator==(const NotDefaultConstructible& other) const
    {
        return other.i == i;
    }
};

struct NotDefaultConstructibleParam {
    using Type = NotDefaultConstructible;
    static Type get_default()
    {
        return NotDefaultConstructible(0);
    }
};

struct DefaultConstructible {
    int i;
    DefaultConstructible()
        : i(0)
    {
    }
    DefaultConstructible(int i_)
        : i(i_)
    {
    }
    ~DefaultConstructible()
    {
    }
    bool operator==(const DefaultConstructible& other) const
    {
        return other.i == i;
    }
};

struct DefaultConstructibleParam {
    using Type = DefaultConstructible;
};

struct NoDefaultMemberFunctionParam {
    using Type = DefaultConstructible;
};

struct MoveOnly {
    MoveOnly()                      = default;
    MoveOnly(const MoveOnly& other) = delete;
    MoveOnly& operator=(const MoveOnly& other) = delete;
    MoveOnly(MoveOnly&& other)                 = default;
    MoveOnly& operator=(MoveOnly&& other) = default;
};

struct MoveOnlyParam {
    using Type = MoveOnly;
};

TEST(TestParameterSet, defaultConstructor)
{
    static_assert(!std::is_default_constructible<epi::ParameterSet<NotDefaultConstructibleParam>>::value,
                  "default constructor missing");
    static_assert(std::is_default_constructible<epi::ParameterSet<DefaultConstructibleParam>>::value,
                  "default constructor missing");
    ASSERT_NO_THROW(epi::ParameterSet<NotDefaultConstructibleParam>(epi::DefaultInit{}));
    ASSERT_NO_THROW(epi::ParameterSet<DefaultConstructibleParam>());
}

TEST(TestParameterSet, defaultInitConstructor)
{
    auto params1 = epi::ParameterSet<IntParam1>(epi::DefaultInit());
    ASSERT_EQ(params1.get<IntParam1>(), 1);

    auto params2 = epi::ParameterSet<IntParam1, DoubleParam>(epi::DefaultInit{});
    ASSERT_EQ(params2.get<IntParam1>(), 1);
    ASSERT_EQ(params2.get<DoubleParam>(), 1.0);

    auto params3 = epi::ParameterSet<IntParam1, DoubleParam, IntParam2>(epi::DefaultInit{});
    ASSERT_EQ(params3.get<IntParam1>(), 1);
    ASSERT_EQ(params3.get<DoubleParam>(), 1.0);
    ASSERT_EQ(params3.get<IntParam2>(), 2);

    auto params4 = epi::ParameterSet<NoDefaultMemberFunctionParam>(epi::DefaultInit{});
    ASSERT_EQ(params4.get<NoDefaultMemberFunctionParam>(), DefaultConstructible());
}

TEST(TestParameterSet, setDefault)
{
    auto params1 = epi::ParameterSet<IntParam1, DoubleParam>();
    ASSERT_EQ(params1.get<IntParam1>(), 0);
    params1.set_default<IntParam1>();
    ASSERT_EQ(params1.get<IntParam1>(), 1);
}

TEST(TestParameterSet, explicitInitConstructors)
{
    static_assert(std::is_constructible<epi::ParameterSet<DoubleParam>, double>::value,
                  "explicit initializing constructor missing");
    static_assert(std::is_constructible<epi::ParameterSet<IntParam1, DoubleParam>, int, double>::value,
                  "explicit initializing constructor missing");
    auto params1 = epi::ParameterSet<IntParam1, DoubleParam>{3, 0.5};
    ASSERT_EQ(params1.get<IntParam1>(), 3);
    ASSERT_EQ(params1.get<DoubleParam>(), 0.5);
}

TEST(TestParameterSet, set)
{
    epi::ParameterSet<IntParam1, DoubleParam, IntParam2> params(epi::DefaultInit{});
    params.set<IntParam1>(3);
    ASSERT_EQ(params.get<IntParam1>(), 3);
}

TEST(TestParameterSet, moveOnly)
{
    static_assert(std::is_default_constructible<epi::ParameterSet<MoveOnlyParam>>::value,
                  "move only parameter not default constructible");
    static_assert(std::is_constructible<epi::ParameterSet<MoveOnlyParam>, MoveOnly&&>::value,
                  "move only parameter not move constructible");
    static_assert(std::is_constructible<epi::ParameterSet<MoveOnlyParam>, epi::DefaultInit>::value,
                  "move only parameter not default initializable");

    auto params                 = epi::ParameterSet<MoveOnlyParam>(MoveOnly());
    params.get<MoveOnlyParam>() = MoveOnly();
    params.set<MoveOnlyParam>(MoveOnly());
}

template <class Mock>
struct MockForeachFuncRef {
    template <class Tag>
    void operator()(typename Tag::Type value, Tag tag)
    {
        mock.invoke(value, tag);
    }
    MockForeachFuncRef(Mock& m)
        : mock(m)
    {
    }
    Mock& mock;
};

TEST(TestParameterSet, foreach)
{
    using testing::An;
    epi::ParameterSet<IntParam1, DoubleParam, IntParam2> params(epi::DefaultInit{});

    struct MockForeachFunc {
        MOCK_METHOD(void, invoke, (int, IntParam1), ());
        MOCK_METHOD(void, invoke, (int, IntParam2), ());
        MOCK_METHOD(void, invoke, (double, DoubleParam), ());
    };

    MockForeachFunc mock;
    {
        testing::InSequence seq{};
        EXPECT_CALL(mock, invoke(An<int>(), An<IntParam1>())).Times(1);
        EXPECT_CALL(mock, invoke(An<double>(), An<DoubleParam>())).Times(1);
        EXPECT_CALL(mock, invoke(An<int>(), An<IntParam2>())).Times(1);
    }

    epi::foreach (params, MockForeachFuncRef<MockForeachFunc>{mock});
}

template <class Mock>
struct MockForeachTagFuncRef {
    template <class Tag>
    void operator()(Tag tag)
    {
        mock.invoke(tag);
    }
    MockForeachTagFuncRef(Mock& m)
        : mock(m)
    {
    }
    Mock& mock;
};

TEST(TestParameterSet, foreach_tag)
{
    epi::ParameterSet<IntParam1, DoubleParam> params;
    struct MockForeachTagFunc {
        MOCK_METHOD(void, invoke, (IntParam1), ());
        MOCK_METHOD(void, invoke, (DoubleParam), ());
    };
    MockForeachTagFunc mock;
    {
        testing::InSequence seq{};
        EXPECT_CALL(mock, invoke(testing::An<IntParam1>())).Times(1);
        EXPECT_CALL(mock, invoke(testing::An<DoubleParam>())).Times(1);
    }
    epi::foreach_tag<epi::ParameterSet<IntParam1, DoubleParam>>(MockForeachTagFuncRef<MockForeachTagFunc>{mock});
}

struct ConstTypeParam {
    using Type = const double;
    static constexpr double get_default()
    {
        return 1.0;
    }
};

TEST(TestParameterSet, constType)
{
    static_assert(std::is_default_constructible<epi::ParameterSet<ConstTypeParam>>::value,
                  "const parameter not default constructible");
    static_assert(std::is_constructible<epi::ParameterSet<ConstTypeParam>, double>::value,
                  "const parameter not constructible");
    static_assert(std::is_constructible<epi::ParameterSet<ConstTypeParam>, epi::DefaultInit>::value,
                  "const parameter not default initializable");
}

TEST(TestParameterSet, equality)
{
    epi::ParameterSet<IntParam1, DoubleParam> a(1, 0.5);
    epi::ParameterSet<IntParam1, DoubleParam> b(1, 0.5);
    epi::ParameterSet<IntParam1, DoubleParam> c(1, 0.6);
    epi::ParameterSet<IntParam1, DoubleParam> d(2, 0.5);
    
    ASSERT_TRUE(a == b);
    ASSERT_TRUE(a != c);
    ASSERT_TRUE(a != d);
}