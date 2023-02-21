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
#include "memilio/geography/regions.h"
#include "gtest/gtest.h"

TEST(TestRegions, get_holidays)
{
    auto a = mio::regions::get_holidays(mio::regions::StateId(9), mio::Date(2020, 10, 15), mio::Date(2020, 11, 15));
    ASSERT_EQ(a.size(), 1);
    ASSERT_EQ(a[0], std::make_pair(mio::Date(2020, 10, 31), mio::Date(2020, 11, 7)));

    auto b = mio::regions::get_holidays(mio::regions::StateId(3), mio::Date(2020, 7, 30), mio::Date(2020, 12, 31));
    ASSERT_EQ(b.size(), 3);
    ASSERT_EQ(b[0], std::make_pair(mio::Date(2020, 7, 16), mio::Date(2020, 8, 27)));
    ASSERT_EQ(b[1], std::make_pair(mio::Date(2020, 10, 12), mio::Date(2020, 10, 24)));
    ASSERT_EQ(b[2], std::make_pair(mio::Date(2020, 12, 23), mio::Date(2021, 1, 9)));
}

TEST(TestRegions, get_state_id)
{
    ASSERT_EQ(mio::regions::get_state_id(int(mio::regions::CountyId(1001))), mio::regions::StateId(1));
    ASSERT_EQ(mio::regions::get_state_id(int(mio::regions::CountyId(2000))), mio::regions::StateId(2));
    ASSERT_EQ(mio::regions::get_state_id(int(mio::regions::CountyId(5970))), mio::regions::StateId(5));
    ASSERT_EQ(mio::regions::get_state_id(int(mio::regions::CountyId(9161))), mio::regions::StateId(9));
}
