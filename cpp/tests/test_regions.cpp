/*
* Copyright (C) 2020-2025 MEmilio
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
#include "memilio/geography/regions.h"
#include "gtest/gtest.h"

TEST(TestRegions, get_holidays)
{

    auto s1 = mio::regions::get_holidays(mio::regions::StateId(1), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s1.size(), 1);
    ASSERT_EQ(s1[0], std::make_pair(mio::Date(2020, 10, 5), mio::Date(2020, 10, 18)));

    auto s2 = mio::regions::get_holidays(mio::regions::StateId(2), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s2.size(), 1);
    ASSERT_EQ(s2[0], std::make_pair(mio::Date(2020, 10, 5), mio::Date(2020, 10, 17)));

    auto s4 = mio::regions::get_holidays(mio::regions::StateId(4), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s4.size(), 1);
    ASSERT_EQ(s4[0], std::make_pair(mio::Date(2020, 10, 12), mio::Date(2020, 10, 25)));

    auto s5 = mio::regions::get_holidays(mio::regions::StateId(5), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s5.size(), 1);
    ASSERT_EQ(s5[0], std::make_pair(mio::Date(2020, 10, 12), mio::Date(2020, 10, 25)));

    auto s6 = mio::regions::get_holidays(mio::regions::StateId(6), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s6.size(), 1);
    ASSERT_EQ(s6[0], std::make_pair(mio::Date(2020, 10, 5), mio::Date(2020, 10, 18)));

    auto s7 = mio::regions::get_holidays(mio::regions::StateId(7), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s7.size(), 1);
    ASSERT_EQ(s7[0], std::make_pair(mio::Date(2020, 10, 12), mio::Date(2020, 10, 24)));

    auto s8 = mio::regions::get_holidays(mio::regions::StateId(8), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s8.size(), 1);
    ASSERT_EQ(s8[0], std::make_pair(mio::Date(2020, 10, 26), mio::Date(2020, 10, 31)));

    auto s9 = mio::regions::get_holidays(mio::regions::StateId(9), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s9.size(), 1);
    ASSERT_EQ(s9[0], std::make_pair(mio::Date(2020, 10, 31), mio::Date(2020, 11, 7)));

    auto s10 = mio::regions::get_holidays(mio::regions::StateId(10), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s10.size(), 1);
    ASSERT_EQ(s10[0], std::make_pair(mio::Date(2020, 10, 12), mio::Date(2020, 10, 24)));

    auto s11 = mio::regions::get_holidays(mio::regions::StateId(11), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s11.size(), 1);
    ASSERT_EQ(s11[0], std::make_pair(mio::Date(2020, 10, 12), mio::Date(2020, 10, 25)));

    auto s12 = mio::regions::get_holidays(mio::regions::StateId(12), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s12.size(), 1);
    ASSERT_EQ(s12[0], std::make_pair(mio::Date(2020, 10, 12), mio::Date(2020, 10, 25)));

    auto s13 = mio::regions::get_holidays(mio::regions::StateId(13), mio::Date(2020, 10, 1), mio::Date(2020, 11, 15));
    ASSERT_EQ(s13.size(), 1);
    ASSERT_EQ(s13[0], std::make_pair(mio::Date(2020, 10, 5), mio::Date(2020, 10, 11)));

    auto s14 = mio::regions::get_holidays(mio::regions::StateId(14), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s14.size(), 1);
    ASSERT_EQ(s14[0], std::make_pair(mio::Date(2020, 10, 19), mio::Date(2020, 11, 1)));

    auto s15 = mio::regions::get_holidays(mio::regions::StateId(15), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s15.size(), 1);
    ASSERT_EQ(s15[0], std::make_pair(mio::Date(2020, 10, 19), mio::Date(2020, 10, 25)));

    auto s16 = mio::regions::get_holidays(mio::regions::StateId(16), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(s16.size(), 1);
    ASSERT_EQ(s16[0], std::make_pair(mio::Date(2020, 10, 17), mio::Date(2020, 10, 31)));

    auto s3 = mio::regions::get_holidays(mio::regions::StateId(3), mio::Date(2020, 7, 30), mio::Date(2020, 12, 31));
    ASSERT_EQ(s3.size(), 3);
    ASSERT_EQ(s3[0], std::make_pair(mio::Date(2020, 7, 16), mio::Date(2020, 8, 27)));
    ASSERT_EQ(s3[1], std::make_pair(mio::Date(2020, 10, 12), mio::Date(2020, 10, 24)));
    ASSERT_EQ(s3[2], std::make_pair(mio::Date(2020, 12, 23), mio::Date(2021, 1, 9)));
}

TEST(TestRegions, get_state_id)
{
    ASSERT_EQ(mio::regions::get_state_id(int(mio::regions::CountyId(1001))), mio::regions::StateId(1));
    ASSERT_EQ(mio::regions::get_state_id(int(mio::regions::CountyId(2000))), mio::regions::StateId(2));
    ASSERT_EQ(mio::regions::get_state_id(int(mio::regions::CountyId(5970))), mio::regions::StateId(5));
    ASSERT_EQ(mio::regions::get_state_id(int(mio::regions::CountyId(9161))), mio::regions::StateId(9));
}
