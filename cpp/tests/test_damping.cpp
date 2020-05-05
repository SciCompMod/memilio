#include "epidemiology/damping.h"
#include <gtest/gtest.h>

TEST(TestDamping, initialDampingIsIdentityEverywhere)
{
    epi::Dampings dampings;
    for (auto x : {-1e100, -12.35, -1e-23, 0.0, 1e-76, 5.67, 1e75}) {
        EXPECT_EQ(dampings.get_factor(0), 1);
    }
}

TEST(TestDamping, dampingContinuesConstantOnBothSides)
{
    double d = 13.4;
    double v = 5.723;

    epi::Dampings dampings;
    dampings.add({d, v});
    for (auto x : {1e-76, 5.67, 1e75}) {
        EXPECT_EQ(dampings.get_factor(0 - x), 1) << "extrapolating before first";
        EXPECT_EQ(dampings.get_factor(d + x), v) << "extrapolating after last";
    }
}

TEST(TestDamping, dampingIsConstantBetweenTwoPoints)
{
    epi::Dampings dampings;
    dampings.add({1, 3.4});
    dampings.add({2, 1.5});

    double eps = 1e-15;
    EXPECT_EQ(dampings.get_factor(1 + eps), 3.4);
    EXPECT_EQ(dampings.get_factor(2 - eps), 3.4);
    EXPECT_EQ(dampings.get_factor(1.5565), 3.4);
}

TEST(TestDamping, dampingJumpsAreClosedOnTheRightSide)
{
    epi::Dampings dampings;
    dampings.add({-1, 3.14});
    dampings.add({2, 2.13});

    const double eps = 1e-15;
    EXPECT_EQ(dampings.get_factor(0), 1);
    EXPECT_EQ(dampings.get_factor(0 - eps), 3.14);
    EXPECT_EQ(dampings.get_factor(2), 2.13);
    EXPECT_EQ(dampings.get_factor(2 - eps), 1);
}

// TEST(TestDamping, duplicatePointsAreUndefined)
// {
//     epi::Dampings dampings;
//     dampings.add({1, 1.3});
//     dampings.add({1, 5.6});
//     //are duplicate points intentionally undefined?
// }