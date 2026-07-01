/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Kilian Volmer
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

#include "memilio/geography/distance.h"
#include <gtest/gtest.h>

TEST(TestDistance, DistanceConstruction)
{
    // Generate same distance in multiple ways:
    mio::geo::Distance distance(10000);
    mio::geo::Distance distance2 = mio::geo::meters(10000);
    mio::geo::Distance distance3 = mio::geo::kilometers(10);
    // Verify that the set Distances are equal.
    EXPECT_TRUE(distance == distance2);
    EXPECT_TRUE(distance == distance3);
    EXPECT_TRUE(distance2 == distance3);
}

TEST(TestDistance, DistanceConversion)
{
    mio::geo::Distance distance = mio::geo::kilometers(5);
    // Verify that the kilometer distance times 1000 equals the meter distance.
    EXPECT_DOUBLE_EQ(distance.kilometers() * 1000, distance.meters());
}

TEST(TestDistance, DistanceComparisons)
{
    mio::geo::Distance smaller_distance = mio::geo::kilometers(5);
    mio::geo::Distance longer_distance  = mio::geo::kilometers(10);
    // Verify that the comparison operators work as expected.
    EXPECT_TRUE(smaller_distance < longer_distance);
    EXPECT_TRUE(longer_distance > smaller_distance);
    EXPECT_TRUE(smaller_distance <= longer_distance);
    EXPECT_TRUE(longer_distance >= smaller_distance);
    EXPECT_TRUE(smaller_distance <= smaller_distance);
    EXPECT_TRUE(longer_distance >= longer_distance);
    EXPECT_TRUE(smaller_distance != longer_distance);
    EXPECT_FALSE(smaller_distance > longer_distance);
    EXPECT_FALSE(longer_distance < smaller_distance);
    EXPECT_FALSE(smaller_distance >= longer_distance);
    EXPECT_FALSE(longer_distance <= smaller_distance);
}

TEST(TestDistance, DistanceArithmetic)
{
    mio::geo::Distance distance1 = mio::geo::meters(5);
    mio::geo::Distance distance2 = mio::geo::meters(10);
    // Verify that addition and subtraction work as expected.
    EXPECT_DOUBLE_EQ((distance1 + distance2).meters(), 15);
    EXPECT_DOUBLE_EQ((distance2 - distance1).meters(), 5);
    EXPECT_DOUBLE_EQ((distance1 += distance2).meters(), 15);
    EXPECT_DOUBLE_EQ((distance1 -= distance2).meters(), 5);
}
