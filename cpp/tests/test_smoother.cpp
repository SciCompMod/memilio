/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Martin J. Kuehn
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
#include "memilio/math/smoother.h"
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

    double val, oldval = lower_step_value;
    for (size_t i = 0; i < vals_left.size(); i++) {
        val = mio::smoother_cosine(vals_left[i], area_smoothed_left, area_smoothed_right, lower_step_value,
                                   higher_step_value);
        EXPECT_EQ(val, lower_step_value) << vals_left[i];

        EXPECT_GE(val, oldval);
        oldval = val;
    }

    for (size_t i = 0; i < vals_left_int.size(); i++) {
        val = mio::smoother_cosine(vals_left_int[i], area_smoothed_left, area_smoothed_right, lower_step_value,
                                   higher_step_value);
        EXPECT_LT(val, lower_step_value + 0.5 * (higher_step_value - lower_step_value)) << " at " << vals_left_int[i];

        EXPECT_GE(val, oldval);
        oldval = val;
    }

    for (size_t i = 0; i < vals_right_int.size(); i++) {
        val = mio::smoother_cosine(vals_right_int[i], area_smoothed_left, area_smoothed_right, lower_step_value,
                                   higher_step_value);
        EXPECT_GT(val, lower_step_value + 0.5 * (higher_step_value - lower_step_value)) << vals_right_int[i];

        EXPECT_GE(val, oldval);
        oldval = val;
    }

    for (size_t i = 0; i < vals_right.size(); i++) {
        val = mio::smoother_cosine(vals_right[i], area_smoothed_left, area_smoothed_right, lower_step_value,
                                   higher_step_value);
        EXPECT_EQ(val, higher_step_value) << vals_right[i];

        EXPECT_GE(val, oldval);
        oldval = val;
    }
}
