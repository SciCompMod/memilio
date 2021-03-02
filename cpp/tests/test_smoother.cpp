#include "epidemiology/math/smoother.h"
#include <gtest/gtest.h>

TEST(TestSmoother, CosineSmootherCheckValues)
{
    double area_smoothed_left  = 1.43;
    double area_smoothed_right = 10 * 1.43;

    double lower_step_value  = 0.05;
    double higher_step_value = 10 * lower_step_value;

    std::vector<double> vals_left{-1.e17, -17., 1.1, 1.42, 1.43};
    std::vector<double> vals_left_int{1.44, 2.55, 2.65, 3.7, 7.864};
    std::vector<double> vals_right_int{7.866, 8.0, 9.8, 12., 14.2};
    std::vector<double> vals_right{14.3, 14.4, 17., 200., 1.e20};

    double val, oldval;
    for (size_t i = 0; i < vals_left.size(); i++) {
        val = epi::smoother_cosine(vals_left[i], area_smoothed_left, area_smoothed_right, lower_step_value,
                                   higher_step_value);
        EXPECT_EQ(val, lower_step_value) << vals_left[i];
        oldval = val;
        EXPECT_GE(val, oldval);
    }

    for (size_t i = 0; i < vals_left_int.size(); i++) {
        val = epi::smoother_cosine(vals_left_int[i], area_smoothed_left, area_smoothed_right, lower_step_value,
                                   higher_step_value);
        EXPECT_LT(val, lower_step_value + 0.5 * (higher_step_value - lower_step_value)) << " at " << vals_left_int[i];
        oldval = val;
        EXPECT_GE(val, oldval);
    }

    for (size_t i = 0; i < vals_right_int.size(); i++) {
        val = epi::smoother_cosine(vals_right_int[i], area_smoothed_left, area_smoothed_right, lower_step_value,
                                   higher_step_value);
        EXPECT_GT(val, lower_step_value + 0.5 * (higher_step_value - lower_step_value)) << vals_right_int[i];
        oldval = val;
        EXPECT_GE(val, oldval);
    }

    for (size_t i = 0; i < vals_right.size(); i++) {
        val = epi::smoother_cosine(vals_right[i], area_smoothed_left, area_smoothed_right, lower_step_value,
                                   higher_step_value);
        EXPECT_EQ(val, higher_step_value) << vals_right[i];
        oldval = val;
        EXPECT_GE(val, oldval);
    }
}
