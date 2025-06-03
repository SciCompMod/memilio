/*
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/utils/date.h"
#include <gtest/gtest.h>

TEST(TestDate, init)
{
    auto date = mio::Date(2021, 3, 12);
    EXPECT_EQ(date.year, 2021);
    EXPECT_EQ(date.month, 3);
    EXPECT_EQ(date.day, 12);

#ifndef NDEBUG
    // check if init fails for non-existing dates
    ASSERT_DEATH(mio::Date(2021, 2, 29), ".*");
    ASSERT_DEATH(mio::Date(2021, 4, 31), ".*");
#endif
}

TEST(TestDate, is_leap_year)
{
    EXPECT_TRUE(mio::is_leap_year(2000));
    EXPECT_TRUE(mio::is_leap_year(2020));
    EXPECT_FALSE(mio::is_leap_year(2021));
    EXPECT_FALSE(mio::is_leap_year(2100));
}

TEST(TestDate, get_month_length)
{
    EXPECT_EQ(mio::get_month_length(mio::Date(2022, 8, 12)), 31);
    EXPECT_EQ(mio::get_month_length(mio::Date(2024, 4, 12)), 30);

    EXPECT_EQ(mio::get_month_length(mio::Date(2022, 2, 12)), 28);
    EXPECT_EQ(mio::get_month_length(mio::Date(2000, 2, 12)), 29);
    EXPECT_EQ(mio::get_month_length(mio::Date(2100, 2, 12)), 28);
}

TEST(TestDate, calculate_partial_sum_of_months)
{
    std::array<int, 12> part_sum = {31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365};
    EXPECT_EQ(mio::calculate_partial_sum_of_months(mio::Date(2022, 8, 12)), part_sum);
    part_sum = {31, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335, 366};
    EXPECT_EQ(mio::calculate_partial_sum_of_months(mio::Date(2024, 4, 12)), part_sum);
}

TEST(TestDate, comparison)
{
    EXPECT_EQ(mio::Date(2021, 3, 12), mio::Date(2021, 3, 12));
    EXPECT_NE(mio::Date(2021, 5, 11), mio::Date(2021, 5, 12));
    EXPECT_NE(mio::Date(2021, 5, 11), mio::Date(2021, 6, 11));
    EXPECT_NE(mio::Date(2021, 5, 11), mio::Date(2022, 5, 11));

    EXPECT_TRUE(mio::Date(2020, 5, 10) < mio::Date(2021, 1, 1));
    EXPECT_TRUE(mio::Date(2020, 5, 10) < mio::Date(2020, 6, 1));
    EXPECT_TRUE(mio::Date(2020, 5, 10) < mio::Date(2020, 5, 11));
    EXPECT_FALSE(mio::Date(2021, 5, 10) < mio::Date(2020, 5, 11));
    EXPECT_FALSE(mio::Date(2020, 5, 10) < mio::Date(2020, 5, 10));
    EXPECT_TRUE(mio::Date(2020, 5, 10) <= mio::Date(2020, 5, 10));
    EXPECT_TRUE(mio::Date(2020, 5, 10) <= mio::Date(2020, 5, 11));
    EXPECT_FALSE(mio::Date(2020, 5, 10) <= mio::Date(2020, 5, 9));

    EXPECT_FALSE(mio::Date(2020, 5, 10) > mio::Date(2021, 1, 1));
    EXPECT_FALSE(mio::Date(2020, 5, 10) > mio::Date(2020, 6, 1));
    EXPECT_FALSE(mio::Date(2020, 5, 10) > mio::Date(2020, 5, 11));
    EXPECT_TRUE(mio::Date(2021, 5, 10) > mio::Date(2020, 5, 11));
    EXPECT_FALSE(mio::Date(2020, 5, 10) > mio::Date(2020, 5, 10));
    EXPECT_TRUE(mio::Date(2020, 5, 10) >= mio::Date(2020, 5, 10));
    EXPECT_FALSE(mio::Date(2020, 5, 10) >= mio::Date(2020, 5, 11));
    EXPECT_TRUE(mio::Date(2020, 5, 10) >= mio::Date(2020, 5, 9));
}

TEST(TestDate, offsetByDays)
{
    EXPECT_EQ(mio::offset_date_by_days({2020, 8, 30}, 3), mio::Date(2020, 9, 2));
    EXPECT_EQ(mio::offset_date_by_days({2024, 8, 24}, 6), mio::Date(2024, 8, 30));
    EXPECT_EQ(mio::offset_date_by_days({2024, 8, 24}, -6), mio::Date(2024, 8, 18));
    EXPECT_EQ(mio::offset_date_by_days({2020, 2, 28}, 1), mio::Date(2020, 2, 29));
    EXPECT_EQ(mio::offset_date_by_days({2021, 2, 28}, 1), mio::Date(2021, 3, 1));
    EXPECT_EQ(mio::offset_date_by_days({2021, 2, 3}, -5), mio::Date(2021, 1, 29));
    EXPECT_EQ(mio::offset_date_by_days({2021, 2, 28}, 0), mio::Date(2021, 2, 28));
    EXPECT_EQ(mio::offset_date_by_days({2020, 12, 31}, 1), mio::Date(2021, 1, 1));
    EXPECT_EQ(mio::offset_date_by_days({2021, 12, 31}, 1), mio::Date(2022, 1, 1));
    EXPECT_EQ(mio::offset_date_by_days({2019, 12, 31}, 367), mio::Date(2021, 1, 1));
    EXPECT_EQ(mio::offset_date_by_days({2021, 12, 31}, 366), mio::Date(2023, 1, 1));
    EXPECT_EQ(mio::offset_date_by_days({2021, 1, 1}, -1), mio::Date(2020, 12, 31));
    EXPECT_EQ(mio::offset_date_by_days({2022, 1, 1}, -1), mio::Date(2021, 12, 31));
    EXPECT_EQ(mio::offset_date_by_days({2021, 1, 1}, -367), mio::Date(2019, 12, 31));
    EXPECT_EQ(mio::offset_date_by_days({2022, 1, 1}, -366), mio::Date(2020, 12, 31));
}

TEST(TestDate, parse)
{
    EXPECT_EQ(mio::parse_date("2020.09.02").value(), mio::Date(2020, 9, 2));
    EXPECT_EQ(mio::parse_date("2021.08.30").value(), mio::Date(2021, 8, 30));
    EXPECT_EQ(mio::parse_date("2020.02.29").value(), mio::Date(2020, 2, 29));
    EXPECT_EQ(mio::parse_date("2021.03.01").value(), mio::Date(2021, 3, 1));
    EXPECT_EQ(mio::parse_date("2021.02.28").value(), mio::Date(2021, 2, 28));
}

TEST(TestDate, getDayInYear)
{
    auto day = mio::get_day_in_year({2020, 1, 21});
    EXPECT_EQ(day, 21);

    day = mio::get_day_in_year({2020, 2, 14});
    EXPECT_EQ(day, 45);

    day = mio::get_day_in_year({2020, 3, 2});
    EXPECT_EQ(day, 62);

    day = mio::get_day_in_year({2020, 12, 27});
    EXPECT_EQ(day, 362);

    day = mio::get_day_in_year({2021, 1, 21});
    EXPECT_EQ(day, 21);

    day = mio::get_day_in_year({2021, 2, 14});
    EXPECT_EQ(day, 45);

    day = mio::get_day_in_year({2021, 3, 2});
    EXPECT_EQ(day, 61);

    day = mio::get_day_in_year({2021, 12, 27});
    EXPECT_EQ(day, 361);
}

TEST(TestDate, getOffset)
{
    auto offset = mio::get_offset_in_days({2020, 8, 30}, {2020, 8, 15});
    EXPECT_EQ(offset, 15);

    offset = mio::get_offset_in_days({2020, 8, 30}, {2020, 8, 31});
    EXPECT_EQ(offset, -1);

    offset = mio::get_offset_in_days({2020, 9, 12}, {2020, 8, 30});
    EXPECT_EQ(offset, 13);

    offset = mio::get_offset_in_days({2020, 7, 25}, {2020, 5, 25});
    EXPECT_EQ(offset, 61);

    offset = mio::get_offset_in_days({2021, 1, 3}, {2020, 12, 31});
    EXPECT_EQ(offset, 3);

    offset = mio::get_offset_in_days({2021, 3, 3}, {2020, 12, 29});
    EXPECT_EQ(offset, 64);

    offset = mio::get_offset_in_days({2021, 11, 30}, {2020, 11, 30});
    EXPECT_EQ(offset, 365);

    offset = mio::get_offset_in_days({2025, 11, 30}, {2020, 11, 30});
    EXPECT_EQ(offset, 5 * 365 + 1);

    offset = mio::get_offset_in_days({2019, 11, 30}, {2020, 11, 30});
    EXPECT_EQ(offset, -366);
}

TEST(TestDate, toIsoString)
{
    EXPECT_EQ(mio::Date(2020, 9, 2).to_iso_string(), "2020-09-02");
    EXPECT_EQ(mio::Date(2021, 8, 30).to_iso_string(), "2021-08-30");
    EXPECT_EQ(mio::Date(2021, 3, 4).to_iso_string(), "2021-03-04");
    EXPECT_EQ(mio::Date(2021, 1, 1).to_iso_string(), "2021-01-01");
    EXPECT_EQ(mio::Date(2020, 2, 29).to_iso_string(), "2020-02-29");
}

TEST(TestDate, streamOutput)
{
    std::ostringstream oss1;
    oss1 << mio::Date(2020, 9, 2);
    EXPECT_EQ(oss1.str(), "2020-09-02");

    std::ostringstream oss2;
    oss2 << mio::Date(2021, 8, 30);
    EXPECT_EQ(oss2.str(), "2021-08-30");

    std::ostringstream oss3;
    oss3 << mio::Date(2021, 3, 4);
    EXPECT_EQ(oss3.str(), "2021-03-04");

    std::ostringstream oss4;
    oss4 << mio::Date(2021, 1, 1);
    EXPECT_EQ(oss4.str(), "2021-01-01");

    std::ostringstream oss5;
    oss5 << mio::Date(2020, 2, 29);
    EXPECT_EQ(oss5.str(), "2020-02-29");
}

TEST(TestDate, formatViaFmt)
{
    EXPECT_EQ(fmt::format("{}", mio::Date(2020, 9, 2)), "2020-09-02");
    EXPECT_EQ(fmt::format("{}", mio::Date(2021, 8, 30)), "2021-08-30");
    EXPECT_EQ(fmt::format("{}", mio::Date(2021, 1, 1)), "2021-01-01");
    EXPECT_EQ(fmt::format("{}", mio::Date(2020, 2, 29)), "2020-02-29");
    EXPECT_EQ(fmt::format("{}", mio::Date(2021, 3, 4)), "2021-03-04");
    EXPECT_EQ(fmt::format("Todays date is: {}", mio::Date(2021, 12, 31)), "Todays date is: 2021-12-31");
}