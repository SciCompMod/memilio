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
#include "memilio/utils/type_safe.h"
#include "gtest/gtest.h"

TEST(TypeSafe, init)
{
    DECL_TYPESAFE(int, TS);
    TS ts(3);
    ASSERT_EQ(ts.get(), 3);
}

TEST(TypeSafe, numericOps)
{
    class TS : public mio::TypeSafe<int, TS>,
               public mio::OperatorAdditionSubtraction<TS>,
               public mio::OperatorScalarMultiplicationDivision<TS, int>
    {
    public:
        using mio::TypeSafe<int, TS>::TypeSafe;
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
    class TS : public mio::TypeSafe<int, TS>, public mio::OperatorComparison<TS>
    {
    public:
        using mio::TypeSafe<int, TS>::TypeSafe;
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
