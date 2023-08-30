/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include "memilio/utils/time_series.h"
#include "memilio/utils/stl_util.h"
#include "matchers.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

template <class T>
using TestTimeSeries = ::testing::Test;

using FloatingPointTypes = ::testing::Types<float, double>;

TYPED_TEST_SUITE(TestTimeSeries, FloatingPointTypes);

TYPED_TEST(TestTimeSeries, createEmpty)
{
    mio::TimeSeries<TypeParam> ts(10);
    ASSERT_EQ(ts.get_num_elements(), 10);
    ASSERT_EQ(ts.get_num_rows(), 11);
    ASSERT_EQ(ts.get_num_time_points(), 0);
    ASSERT_EQ(ts.get_capacity(), 0);
}

TYPED_TEST(TestTimeSeries, createInit)
{
    auto v = mio::TimeSeries<TypeParam>::Vector::Random(5).eval();
    mio::TimeSeries<TypeParam> ts(0.0, v);
    ASSERT_EQ(ts.get_num_elements(), 5);
    ASSERT_EQ(ts.get_num_rows(), 6);
    ASSERT_EQ(ts.get_num_time_points(), 1);
    ASSERT_EQ(ts.get_capacity(), 1);
    ASSERT_EQ(ts.get_time(0), 0.0);
    ASSERT_EQ(print_wrap(ts.get_value(0)), print_wrap(v));
}

TYPED_TEST(TestTimeSeries, zeroElements)
{
    mio::TimeSeries<TypeParam> ts(0);
    ASSERT_EQ(ts.get_num_elements(), 0);
    ASSERT_EQ(ts.get_num_rows(), 1);
    ASSERT_EQ(ts.get_num_time_points(), 0);
    ASSERT_EQ(ts.get_capacity(), 0);
    ASSERT_NO_FATAL_FAILURE(ts.add_time_point(0.0));
}

TYPED_TEST(TestTimeSeries, addPoints)
{
    mio::TimeSeries<TypeParam> ts(5);
    ts.add_time_point(0.0);
    ASSERT_EQ(ts.get_num_time_points(), 1);
    ASSERT_EQ(ts.get_capacity(), 1 << 0);
    ts.add_time_point(1.0);
    ASSERT_EQ(ts.get_num_time_points(), 2);
    ASSERT_EQ(ts.get_capacity(), 1 << 1);
    ts.add_time_point(2.0);
    ASSERT_EQ(ts.get_num_time_points(), 3);
    ASSERT_EQ(ts.get_capacity(), 1 << 2);

    float i = 3.0f;
    while (i < 7) {
        ts.add_time_point(i);
        ++i;
    }
    ASSERT_EQ(ts.get_num_time_points(), 7);
    ASSERT_EQ(ts.get_capacity(), 1 << 3);

    while (i < 123456) {
        ts.add_time_point(i);
        ++i;
    }
    ASSERT_EQ(ts.get_num_time_points(), 123456);
    ASSERT_EQ(ts.get_capacity(), 1 << 17);
}

TYPED_TEST(TestTimeSeries, assignValues)
{
    mio::TimeSeries<TypeParam> ts(2);
    auto v0               = mio::TimeSeries<TypeParam>::Vector::Random(2).eval();
    ts.add_time_point(0.) = v0;
    auto v1               = mio::TimeSeries<TypeParam>::Vector::Random(2).eval();
    ts.add_time_point(1.);
    ts[1] = v1;
    ts.add_time_point(2., mio::TimeSeries<TypeParam>::Vector::Constant(2, 1));
    auto v2 = mio::TimeSeries<TypeParam>::Vector::Constant(2, 1);

    ASSERT_EQ(print_wrap(ts[0]), print_wrap(v0));
    ASSERT_EQ(print_wrap(ts[1]), print_wrap(v1));
    ASSERT_EQ(print_wrap(ts[2]), print_wrap(v2));
}

TYPED_TEST(TestTimeSeries, copyEmpty)
{
    mio::TimeSeries<TypeParam> ts(10);
    mio::TimeSeries<TypeParam> ts_copy1(ts);
    mio::TimeSeries<TypeParam> ts_copy2(1);
    ts_copy2 = ts;

    for (auto&& ts_copy : {&ts_copy1, &ts_copy2}) {
        ASSERT_EQ(ts_copy->get_num_elements(), 10);
        ASSERT_EQ(ts_copy->get_num_rows(), 11);
        ASSERT_EQ(ts_copy->get_num_time_points(), 0);
        ASSERT_EQ(ts_copy->get_capacity(), 0);
    }
}

TYPED_TEST(TestTimeSeries, reserve)
{
    mio::TimeSeries<TypeParam> ts(2);
    ts.reserve(10);
    ASSERT_EQ(ts.get_capacity(), 16);
    ts.reserve(200);
    ASSERT_EQ(ts.get_capacity(), 256);
    ts.reserve(10);
    ASSERT_EQ(ts.get_capacity(), 256);
}

TYPED_TEST(TestTimeSeries, constAccess)
{
    mio::TimeSeries<TypeParam> ts(1);
    ts.add_time_point(0., mio::TimeSeries<TypeParam>::Vector::Random(1));
    const auto& constref = ts;
    static_assert(
        std::is_same<decltype(constref[0]), Eigen::Ref<const typename mio::TimeSeries<TypeParam>::Vector>>::value,
        "wrong type");
    static_assert(std::is_same<decltype(constref.get_value(0)),
                               Eigen::Ref<const typename mio::TimeSeries<TypeParam>::Vector>>::value,
                  "wrong type");
    static_assert(std::is_same<decltype(constref.get_last_value()),
                               Eigen::Ref<const typename mio::TimeSeries<TypeParam>::Vector>>::value,
                  "wrong type");
    ASSERT_EQ(print_wrap(ts[0]), print_wrap(constref[0]));
}

#ifndef NDEBUG
TYPED_TEST(TestTimeSeries, createInvalidDim)
{
    if (std::is_signed<Eigen::Index>::value) {
        auto create = []() {
            mio::TimeSeries<TypeParam> ts(-1);
            return ts;
        };
        ASSERT_DEATH(create(), ".*");
    }
}

TYPED_TEST(TestTimeSeries, accessInvalidRange)
{
    mio::TimeSeries<TypeParam> ts(1);
    for (Eigen::Index i = 0; i < 123; i++) {
        ts.add_time_point();
    }
    ASSERT_DEATH(ts.get_value(-1), testing::ContainsRegex(".*"));
    ASSERT_DEATH(ts.get_value(123), testing::ContainsRegex(".*"));
    ASSERT_DEATH(ts.get_value(1231556), testing::ContainsRegex(".*"));
}
#endif

TYPED_TEST(TestTimeSeries, data)
{
    mio::TimeSeries<TypeParam> ts(1);
    auto v0 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 0.5);
    auto v1 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 1.5);
    auto v2 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 2.5);
    auto v3 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 3.5);
    ts.add_time_point(0.0, v0);
    ts.add_time_point(1.0, v1);
    ts.add_time_point(2.0, v2);
    ts.add_time_point(3.0, v3);

    auto data_range = mio::make_range(ts.data(), ts.data() + 8);
    ASSERT_THAT(data_range, testing::ElementsAre(TypeParam(0.0), TypeParam(0.5), TypeParam(1.0), TypeParam(1.5),
                                                 TypeParam(2.0), TypeParam(2.5), TypeParam(3.0), TypeParam(3.5)));
}

TYPED_TEST(TestTimeSeries, iteratorsRange)
{
    mio::TimeSeries<TypeParam> ts(1);
    auto v0 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 0.5);
    auto v1 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 1.5);
    auto v2 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 2.5);
    auto v3 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 3.5);
    ts.add_time_point(0.0, v0);
    ts.add_time_point(1.0, v1);
    ts.add_time_point(2.0, v2);
    ts.add_time_point(3.0, v3);

    //the range-loops and range-assert check the same condition in different ways
    int i = 0;
    for (auto&& v : ts) {
        ASSERT_EQ(print_wrap(v), print_wrap(ts[i]));
        ++i;
    }
    i                       = 0;
    const auto& ts_constref = ts;
    for (auto&& v : ts_constref) {
        ASSERT_EQ(print_wrap(v), print_wrap(ts[i]));
        ++i;
    }
    i = 3;
    for (auto&& v : mio::make_range(ts.rbegin(), ts.rend())) {
        ASSERT_EQ(print_wrap(v), print_wrap(ts[i]));
        --i;
    }
    ASSERT_THAT(ts, testing::ElementsAre(v0, v1, v2, v3));
    ASSERT_THAT(ts_constref, testing::ElementsAre(v0, v1, v2, v3));
    ASSERT_THAT(mio::make_range(ts.rbegin(), ts.rend()), testing::ElementsAre(v3, v2, v1, v0));
}

TYPED_TEST(TestTimeSeries, timeIteratorsRange)
{
    mio::TimeSeries<TypeParam> ts(1);
    auto v0 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 0.5);
    auto v1 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 1.5);
    auto v2 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 2.5);
    auto v3 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 3.5);
    ts.add_time_point(0.0, v0);
    ts.add_time_point(1.0, v1);
    ts.add_time_point(2.0, v2);
    ts.add_time_point(3.0, v3);

    //the range-loops and range-assert check the same condition in different ways
    int i = 0;
    for (auto&& t : ts.get_times()) {
        ASSERT_EQ(t, ts.get_time(i));
        ++i;
    }
    i                       = 0;
    const auto& ts_constref = ts;
    for (auto&& t : ts_constref.get_times()) {
        ASSERT_EQ(t, ts.get_time(i));
        ++i;
    }
    i = 3;
    for (auto&& t : ts.get_reverse_times()) {
        ASSERT_EQ(t, ts.get_time(i));
        --i;
    }
    ASSERT_THAT(ts.get_times(), testing::ElementsAre(TypeParam(0.0), TypeParam(1.0), TypeParam(2.0), TypeParam(3.0)));
    ASSERT_THAT(ts_constref.get_times(),
                testing::ElementsAre(TypeParam(0.0), TypeParam(1.0), TypeParam(2.0), TypeParam(3.0)));
    ASSERT_THAT(ts.get_reverse_times(),
                testing::ElementsAre(TypeParam(3.0), TypeParam(2.0), TypeParam(1.0), TypeParam(0.0)));
}

TYPED_TEST(TestTimeSeries, iteratorsRandomAccess)
{
    mio::TimeSeries<TypeParam> ts(1);
    auto v0 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 0.5);
    auto v1 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 1.5);
    auto v2 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 2.5);
    auto v3 = mio::TimeSeries<TypeParam>::Vector::Constant(1, 3.5);
    ts.add_time_point(0.0, v0);
    ts.add_time_point(1.0, v1);
    ts.add_time_point(2.0, v2);
    ts.add_time_point(3.0, v3);

    auto it0   = ts.begin();
    auto it1   = ts.begin() + 1;
    auto it2   = ts.begin() + 2;
    auto it3   = ts.begin() + 3;
    auto itEnd = ts.end();

    //deref
    ASSERT_EQ(print_wrap(*it0), print_wrap(v0));
    ASSERT_EQ(print_wrap(*it1), print_wrap(v1));
    ASSERT_EQ(print_wrap(*it2), print_wrap(v2));
    ASSERT_EQ(print_wrap(*it3), print_wrap(v3));

    //index
    ASSERT_EQ(print_wrap(it1[1]), print_wrap(v2));

    //addition
    auto it2c1 = it2;
    ASSERT_EQ(++it2c1, it3);
    auto it2c2 = it2;
    ASSERT_EQ(it2c2++, it2);
    auto it2c3 = it2;
    ASSERT_EQ(it2c3 += 1, it3);
    ASSERT_EQ(it2 + 2, itEnd);

    //subtraction
    auto it3c1 = it3;
    ASSERT_EQ(--it3c1, it2);
    auto it3c2 = it3;
    ASSERT_EQ(it3c2--, it3);
    auto it3c3 = it3;
    ASSERT_EQ(it3c3 -= 3, it0);
    ASSERT_EQ(itEnd - 1, it3);

    //comparison
    ASSERT_LT(it0, it1);
    ASSERT_GT(it3, it1);
    ASSERT_LE(it1, itEnd);
    ASSERT_GE(it2, it0);
}

TYPED_TEST(TestTimeSeries, create)
{
    auto ts = mio::TimeSeries<double>::zero(5, 10);
    for (int i = 0; i < 5; i++) {
        ASSERT_EQ(ts.get_time(i), 0.0);
        for (int j = 0; j < 10; j++) {
            ASSERT_EQ(ts[i][j], 0.0);
        }
    }
}
