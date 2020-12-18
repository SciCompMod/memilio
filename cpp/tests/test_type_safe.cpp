#include "epidemiology/utils/type_safe.h"
#include "gtest/gtest.h"

TEST(TypeSafe, init)
{
    DECL_TYPESAFE(TS, int);
    TS ts(3);
    ASSERT_EQ(ts.get(), 3);
}

TEST(TypeSafe, numericOps)
{
    class TS : public epi::TypeSafe<int, TS>,
               public epi::OperatorAddSub<TS>,
               public epi::OperatorScalarMultDiv<TS, int>,
               public epi::OperatorIncrDecr<TS>
    {
    public:
        using epi::TypeSafe<int, TS>::TypeSafe;
    };

    {
        TS ts1(3), ts2(2);
        ASSERT_EQ((ts1 + ts2).get(), 5);
        ASSERT_EQ((ts1 += ts2).get(), 5);
    }

    {
        TS ts1(3), ts2(2);
        ASSERT_EQ((ts1 - ts2).get(), 1);
        ASSERT_EQ((ts1 -= ts2).get(), 1);
    }

    {
        TS ts1(3);
        ASSERT_EQ((ts1 * 2).get(), 6);
        ASSERT_EQ((ts1 *= 2).get(), 6);
    }
    {
        TS ts1(3);
        ASSERT_EQ((ts1 / 2).get(), 1);
        ASSERT_EQ((ts1 /= 2).get(), 1);
    }

    {
        TS ts1(3);
        ASSERT_EQ((++ts1).get(), 4);
        ASSERT_EQ((ts1++).get(), 4);
        ASSERT_EQ(ts1.get(), 5);
    }

    {
        TS ts1(3);
        ASSERT_EQ((--ts1).get(), 2);
        ASSERT_EQ((ts1--).get(), 2);
        ASSERT_EQ(ts1.get(), 1);
    }
}

TEST(TypeSafe, comparisonOps)
{
    class TS : public epi::TypeSafe<int, TS>, public epi::OperatorComparison<TS>
    {
    public:
        using epi::TypeSafe<int, TS>::TypeSafe;
    };

    TS ts1(3), ts2(2), ts3(3);

    ASSERT_NE(ts1, ts2);
    ASSERT_EQ(ts1, ts3);
    ASSERT_LT(ts2, ts1);
    ASSERT_LE(ts2, ts1);
    ASSERT_LE(ts3, ts1);
    ASSERT_GT(ts1, ts2);
    ASSERT_GE(ts1, ts2);
    ASSERT_GE(ts1, ts3);
}