/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Jan Kleinert
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
#include "memilio/utils/logging.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/custom_index_array.h"
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
    static_assert(!std::is_constructible<mio::ParameterSet<NotDefaultConstructibleParam>, mio::NoDefaultInit>::value,
                  "default constructor missing");
    static_assert(std::is_default_constructible<mio::ParameterSet<DefaultConstructibleParam>>::value,
                  "default constructor missing");
    ASSERT_NO_THROW(mio::ParameterSet<NotDefaultConstructibleParam>());
    ASSERT_NO_THROW(mio::ParameterSet<DefaultConstructibleParam>(mio::NoDefaultInit{}));
}

TEST(TestParameterSet, defaultInitConstructor)
{
    auto params1 = mio::ParameterSet<IntParam1>();
    ASSERT_EQ(params1.get<IntParam1>(), 1);

    auto params2 = mio::ParameterSet<IntParam1, DoubleParam>();
    ASSERT_EQ(params2.get<IntParam1>(), 1);
    ASSERT_EQ(params2.get<DoubleParam>(), 1.0);

    auto params3 = mio::ParameterSet<IntParam1, DoubleParam, IntParam2>();
    ASSERT_EQ(params3.get<IntParam1>(), 1);
    ASSERT_EQ(params3.get<DoubleParam>(), 1.0);
    ASSERT_EQ(params3.get<IntParam2>(), 2);

    auto params4 = mio::ParameterSet<NoDefaultMemberFunctionParam>();
    ASSERT_EQ(params4.get<NoDefaultMemberFunctionParam>(), DefaultConstructible());
}

TEST(TestParameterSet, setDefault)
{
    auto params1 = mio::ParameterSet<IntParam1, DoubleParam>(mio::NoDefaultInit{});
    ASSERT_EQ(params1.get<IntParam1>(), 0);
    params1.set_default<IntParam1>();
    ASSERT_EQ(params1.get<IntParam1>(), 1);
}

TEST(TestParameterSet, set)
{
    mio::ParameterSet<IntParam1, DoubleParam, IntParam2> params;
    params.set<IntParam1>(3);
    ASSERT_EQ(params.get<IntParam1>(), 3);
}

TEST(TestParameterSet, moveOnly)
{
    static_assert(std::is_default_constructible<mio::ParameterSet<MoveOnlyParam>>::value,
                  "move only parameter not default constructible");
    static_assert(std::is_constructible<mio::ParameterSet<MoveOnlyParam>>::value,
                  "move only parameter not default initializable");

    mio::ParameterSet<MoveOnlyParam> params;
    params.set<MoveOnlyParam>(MoveOnly());
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

TEST(TestParameterSet, customIndexArray)
{
    struct AgeGroup {
    };
    struct Income {
    };

    struct ParamType1 {
        using Type = mio::CustomIndexArray<double, AgeGroup>;
        static Type get_default(mio::Index<AgeGroup> n_agegroups, mio::Index<Income>)
        {
            return Type({n_agegroups}, 0.5);
        }
    };

    struct ParamType2 {
        using Type = mio::CustomIndexArray<int, AgeGroup, Income>;
        static Type get_default(mio::Index<AgeGroup> n_agegroups, mio::Index<Income> n_incomegroups)
        {
            return Type({n_agegroups, n_incomegroups}, 42);
        }
    };

    auto params = mio::ParameterSet<ParamType1, ParamType2>(mio::Index<AgeGroup>(2), mio::Index<Income>(3));
    params.get<ParamType1>()[{mio::Index<AgeGroup>(0)}] = 0.5;
    params.get<ParamType1>()[{mio::Index<AgeGroup>(1)}] = 1.5;
    EXPECT_NEAR(params.get<ParamType1>()[{mio::Index<AgeGroup>(0)}], 0.5, 1e-14);
    EXPECT_NEAR(params.get<ParamType1>()[{mio::Index<AgeGroup>(1)}], 1.5, 1e-14);
    EXPECT_EQ((params.get<ParamType2>()[{mio::Index<AgeGroup>(0), mio::Index<Income>(0)}]), 42);
    EXPECT_EQ((params.get<ParamType2>()[{mio::Index<AgeGroup>(0), mio::Index<Income>(1)}]), 42);
    params.get<ParamType2>()[{mio::Index<AgeGroup>(0), mio::Index<Income>(1)}] = -42;
    EXPECT_EQ((params.get<ParamType2>()[{mio::Index<AgeGroup>(0), mio::Index<Income>(1)}]), -42);
}

TEST(TestParameterSet, foreach)
{
    using testing::An;
    mio::ParameterSet<IntParam1, DoubleParam, IntParam2> params;

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

    mio::foreach (params, MockForeachFuncRef<MockForeachFunc>{mock});
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
    mio::ParameterSet<IntParam1, DoubleParam> params;
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
    mio::foreach_tag<mio::ParameterSet<IntParam1, DoubleParam>>(MockForeachTagFuncRef<MockForeachTagFunc>{mock});
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
    static_assert(std::is_default_constructible<mio::ParameterSet<ConstTypeParam>>::value,
                  "const parameter not default constructible");
    static_assert(std::is_constructible<mio::ParameterSet<ConstTypeParam>>::value,
                  "const parameter not default initializable");
}

TEST(TestParameterSet, equality)
{
    mio::ParameterSet<IntParam1, DoubleParam> a;
    a.set<IntParam1>(1);
    a.set<DoubleParam>(0.5);

    mio::ParameterSet<IntParam1, DoubleParam> b;
    b.set<IntParam1>(1);
    b.set<DoubleParam>(0.5);

    mio::ParameterSet<IntParam1, DoubleParam> c;
    c.set<IntParam1>(1);
    c.set<DoubleParam>(0.6);

    mio::ParameterSet<IntParam1, DoubleParam> d;
    d.set<IntParam1>(2);
    d.set<DoubleParam>(0.5);

    ASSERT_TRUE(a == b);
    ASSERT_TRUE(a != c);
    ASSERT_TRUE(a != d);
}
