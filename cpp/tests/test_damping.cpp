#include "epidemiology/secir/damping.h"
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

    for (auto x : {1.00001, 5.67, 1e75}) {
        EXPECT_EQ(dampings.get_factor(0 - x), 1) << "extrapolating before first";
        EXPECT_EQ(dampings.get_factor(d + x), v) << "extrapolating after last";
    }
}

TEST(TestDamping, dampingIsConstantBetweenTwoPoints)
{
    epi::Dampings dampings;
    dampings.add({1, 3.4});
    dampings.add({12, 1.5});

    double eps = 1e-15;
    EXPECT_EQ(dampings.get_factor(1 + 1 + eps), 3.4);
    EXPECT_EQ(dampings.get_factor(12 - 1 - eps), 3.4);
    EXPECT_EQ(dampings.get_factor(6.5565), 3.4);
}

TEST(TestDamping, dampingJumpsMonotonouslySmoothed)
{
    epi::Dampings dampings;
    dampings.add({2, 3.14});
    dampings.add({3, 2.13});
    dampings.add({3.5, 0.13});
    dampings.add({4.5, 4.13});

    int nb_steps = 20;
    double step  = 1.0 / (double)nb_steps;
    for (int i = 0; i < nb_steps; i++) {
        EXPECT_EQ(dampings.get_factor(2 + i * step) > dampings.get_factor(2 + (i + 1) * step), true);
        // printf("\n %f %f ", d - 0.5 + (double)i / 100, dampings.get_factor(d - 0.5 + (double)i / 100));
    }

    for (int i = 0; i < nb_steps; i++) {
        EXPECT_EQ(dampings.get_factor(3 + 0.5 * i * step) > dampings.get_factor(3 + 0.5 * (i + 1) * step), true);
    }

    for (int i = 0; i < nb_steps; i++) {
        EXPECT_EQ(dampings.get_factor(3.5 + i * step) < dampings.get_factor(3.5 + (i + 1) * step), true);
    }
}

TEST(TestDamping, duplicatePointsOverwriteTheOldPoint)
{
    epi::Dampings dampings;
    dampings.add({2, 1.3});
    dampings.add({2, 0.01});
    dampings.add({2, 5.6});
    dampings.add({5, 2.5});

    //old value is overwritten
    EXPECT_EQ(dampings.get_factor(2), 5.6);
    EXPECT_EQ(dampings.get_factor(3), 5.6);

    //other values are unaffected
    EXPECT_EQ(dampings.get_factor(0.5), 1);
    EXPECT_EQ(dampings.get_factor(5), 2.5);
    EXPECT_EQ(dampings.get_factor(5.0341), 2.5);
}