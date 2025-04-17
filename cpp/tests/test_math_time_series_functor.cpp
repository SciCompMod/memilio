/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: René Schmieding
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
#include "memilio/math/time_series_functor.h"
#include "random_number_test.h"

#include "gtest/gtest.h"
#include <vector>

using TestMathTimeSeriesFunctor = RandomNumberTest;

TEST_F(TestMathTimeSeriesFunctor, zero)
{
    // Test that the default initialized functor always returns zero, using a random evaluation point.

    const int num_evals = 100;

    // initialize functor using the default ctor
    mio::TimeSeriesFunctor<double> tsf;

    // check one deterministic value first to avoid flooding the test output with failed tests
    ASSERT_EQ(tsf(0.0), 0.0);
    // verify output
    for (int i = 0; i < num_evals; i++) {
        auto random_t_eval = this->random_number();
        EXPECT_EQ(tsf(random_t_eval), 0.0);
    }
}

TEST_F(TestMathTimeSeriesFunctor, linearInterpolation)
{
    // Test that linear interpolation works for a piecewise linear function.

    // continuous function that is constant 1 for t<0, linear in [0, 1] with slope 2, and constant 3 for t>1
    const auto pcw_lin_fct = [&](double t) {
        return 1 + 2 * std::clamp(t, 0.0, 1.0); // this looks like .../```
    };

    mio::TimeSeriesFunctor<double> tsf(mio::TimeSeriesFunctorType::LinearInterpolation, {{0., 1.}, {1., 3.}});

    // go from -1/4 to 5/4 in steps of size 1/4, with values 1.0, 1.0, 1.5, 2.0, 2.5, 3.0, 3.0
    for (double t = -0.25; t < 1.3; t += 0.25) {
        EXPECT_NEAR(tsf(t), pcw_lin_fct(t), 1e-14);
    }
}

TEST_F(TestMathTimeSeriesFunctor, linearInterpolationRandomized)
{
    // Test that the LinearInterpolation-functor correctly reproduces a piecewise linear function, using random
    // samples. Since the initialization uses unsorted data, this also checks that the data gets sorted
    const int num_evals = 1000;

    const double t_min = -1, t_max = 1, t_mid = this->random_number(t_min, t_max);
    const double slope1 = this->random_number();
    const double slope2 = this->random_number();
    const double height = this->random_number();

    // continuous function with different slopes between t_min, t_mid and t_max, constant otherwise
    const auto pcw_lin_fct = [&](double t) {
        return height + slope1 * std::clamp(t - t_min, 0.0, t_mid - t_min) +
               slope2 * std::clamp(t - t_mid, 0.0, t_max - t_mid);
    };

    // initialize the data with the critical points
    std::vector<std::vector<double>> unsorted_data{
        {t_max, pcw_lin_fct(t_max)}, {t_min, pcw_lin_fct(t_min)}, {t_mid, pcw_lin_fct(t_mid)}};
    // randomly add a few more evaluations in between
    for (int i = 0; i < 10; i++) {
        const double t = this->random_number(-1.0, 1.0);
        unsorted_data.push_back({t, pcw_lin_fct(t)});
    }
    // initialize functor
    mio::TimeSeriesFunctor<double> tsf(mio::TimeSeriesFunctorType::LinearInterpolation, unsorted_data);

    // check one deterministic value first to avoid flooding the test output with failed tests
    ASSERT_NEAR(tsf(0.5 * (t_max - t_min)), pcw_lin_fct(0.5 * (t_max - t_min)), 1e-10);
    // verify output
    for (int i = 0; i < num_evals; i++) {
        // sample in the interval [t_min - (t_max - t_min) / 4, t_max + (t_max - tmin) / 4]
        double random_t_eval = this->random_number(1.25 * t_min - 0.25 * t_max, 1.25 * t_max - 0.25 * t_min);
        EXPECT_NEAR(tsf(random_t_eval), pcw_lin_fct(random_t_eval), 1e-10) << "i = " << i;
    }
}

TEST_F(TestMathTimeSeriesFunctor, unhandledTypes)
{
    GTEST_FLAG_SET(death_test_style, "threadsafe");
    // check that the functor does not accept unhandled types.

    const auto unhandled_type = (mio::TimeSeriesFunctorType)-1;

    // check constructor assert
    EXPECT_DEBUG_DEATH(mio::TimeSeriesFunctor<double>(unhandled_type, mio::TimeSeries<double>(0)),
                       "Unhandled TimeSeriesFunctorType!");

    // abuse default_serialize to set an invalid type
    mio::TimeSeriesFunctor<double> functor;
    std::get<0>(functor.default_serialize().named_refs).value = unhandled_type;

    // check assert in functor call
    EXPECT_DEBUG_DEATH(functor(0.0), "Unhandled TimeSeriesFunctorType!");
}
