/* 
* Copyright (C) 2020-2025 MEmilio
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

#include "memilio/geography/rtree.h"
#include "memilio/geography/geolocation.h"
#include "random_number_test.h"
#include <gtest/gtest.h>

TEST(TestGeography, compareGeographicalLocation)
{
    // Generate GeographicalLocations
    mio::geo::GeographicalLocation geographical_location  = {10.5, 52.2};
    mio::geo::GeographicalLocation geographical_location2 = {10.5, 52.2};
    mio::geo::GeographicalLocation geographical_location3 = {9.5, 52};
    // Verify that the set GeographicalLocations are equal.
    EXPECT_TRUE(geographical_location == geographical_location2);
    // Verify that the set GeographicalLocations are not equal.
    EXPECT_FALSE(geographical_location == geographical_location3);
    // Verify that the set GeographicalLocations are not equal.
    EXPECT_TRUE(geographical_location != geographical_location3);
    // Verify that the set GeographicalLocations are equal.
    EXPECT_FALSE(geographical_location != geographical_location2);
    // Verify that the set GeographicalLocations are close.
    EXPECT_TRUE(geographical_location.is_close(geographical_location2));
    // Verify that the set GeographicalLocations are not close.
    EXPECT_FALSE(geographical_location.is_close(geographical_location3));
    // Verify that the set GeographicalLocations are close with large tolerance.
    EXPECT_TRUE(geographical_location.is_close(geographical_location3, mio::geo::kilometers(600)));
}

TEST(TestGeography, Distance)
{
    // Generate GeographicalLocations
    auto bonn   = mio::geo::GeographicalLocation(50.7, 7.1);
    auto berlin = mio::geo::GeographicalLocation(52.5, 13.4);
    // use negative coordinates
    auto neumayer = mio::geo::GeographicalLocation(-70.6, -8.2);
    // Test that the distance function is symmetric
    EXPECT_DOUBLE_EQ(bonn.distance(berlin).meters(), berlin.distance(bonn).meters());
    mio::geo::Distance distance = mio::geo::kilometers(478.7);
    // case: distance bonn - berlin; expect: ~478.7km
    EXPECT_LT(abs((bonn.distance(berlin) - distance).kilometers()), 0.1);
    // case: distance bonn - bonn; expect: 0km
    EXPECT_DOUBLE_EQ(bonn.distance(bonn).meters(), 0.);
    // case: distance neumayer - bonn; expect: ~13543.7km
    EXPECT_LT(abs(neumayer.distance(bonn).kilometers() - 13543.7), 1.);
}

TEST(TestGeography, rtreeConstructionVector)
{
    // Generate a vector of GeographicalLocations
    std::vector<mio::geo::GeographicalLocation> locations;
    locations.push_back(mio::geo::GeographicalLocation(50.7, 6.0));
    locations.push_back(mio::geo::GeographicalLocation(52.0, 7.0));
    locations.push_back(mio::geo::GeographicalLocation(53.6, 10.2));
    // Generate a RTree object using the vector
    mio::geo::RTree tree(locations);
    // Verify that the size of the tree is equal to the size of the vector
    EXPECT_EQ(tree.size(), locations.size());
}

TEST(TestGeography, rtreeConstructionRange)
{
    // Generate a vector of GeographicalLocations
    std::vector<mio::geo::GeographicalLocation> locations;
    locations.push_back(mio::geo::GeographicalLocation(50.7, 6.0));
    locations.push_back(mio::geo::GeographicalLocation(52.0, 7.0));
    locations.push_back(mio::geo::GeographicalLocation(53.6, 10.3));
    // Generate a RTree object using begin and end iterators
    auto tree = mio::geo::RTree(locations.begin(), locations.end());
    // Verify that the size of the tree is equal to the size of the vector
    EXPECT_EQ(tree.size(), locations.size());
}

TEST(TestGeography, rtreesize)
{
    // Generate a vector of GeographicalLocations
    std::vector<mio::geo::GeographicalLocation> locations;
    locations.push_back(mio::geo::GeographicalLocation(50.7, 6.0));
    locations.push_back(mio::geo::GeographicalLocation(52.0, 7.0));
    locations.push_back(mio::geo::GeographicalLocation(53.6, 10.3));
    // Generate a RTree object
    auto rtree = mio::geo::RTree(locations.begin(), locations.end());
    // Verify that the size of the rtree is 3
    EXPECT_EQ(rtree.size(), 3);
}

TEST(TestGeography, rtreeNN)
{
    // Generate a vector of GeographicalLocations
    std::vector<mio::geo::GeographicalLocation> locations;
    locations.push_back(mio::geo::GeographicalLocation(50.7, 6.0));
    locations.push_back(mio::geo::GeographicalLocation(52.0, 7.0));
    locations.push_back(mio::geo::GeographicalLocation(53.6, 10.2));
    // Generate a RTree object
    auto rtree = mio::geo::RTree(locations.begin(), locations.end());
    // Verify that the nearest neighbor is the first location
    EXPECT_EQ(rtree.nearest_neighbor_indices(mio::geo::GeographicalLocation(50.7, 6.0), 1)[0], 0);
}

TEST(TestGeography, rtreein_range_approx)
{
    // Generate a vector of GeographicalLocations
    std::vector<mio::geo::GeographicalLocation> locations;
    locations.push_back(mio::geo::GeographicalLocation(50.7, 6.0));
    locations.push_back(mio::geo::GeographicalLocation(52.0, 7.0));
    locations.push_back(mio::geo::GeographicalLocation(53.6, 10.2));
    // Generate a RTree object
    auto rtree = mio::geo::RTree(locations.begin(), locations.end());
    // Verify that the in-range queries returns the correct number of results
    EXPECT_EQ(
        rtree.in_range_indices_approximate(mio::geo::GeographicalLocation(50.9, 6.8), mio::geo::kilometers(150)).size(),
        2);
    EXPECT_EQ(rtree.in_range_indices(mio::geo::GeographicalLocation(50.9, 6.8), mio::geo::kilometers(150)).size(), 2);
}

TEST(TestGeography, rtreein_range_multiple_radii)
{
    // Generate a vector of GeographicalLocations
    std::vector<mio::geo::GeographicalLocation> locations;
    locations.push_back(mio::geo::GeographicalLocation(50.7, 6.0));
    locations.push_back(mio::geo::GeographicalLocation(52.0, 7.0));
    locations.push_back(mio::geo::GeographicalLocation(53.6, 10.2));
    // Generate a RTree object
    auto rtree = mio::geo::RTree(locations.begin(), locations.end());
    // Run in_range queries for three different ranges
    auto result =
        rtree.in_range_indices_query(mio::geo::GeographicalLocation(51.4, 7.4),
                                     {mio::geo::kilometers(130), mio::geo::kilometers(320), mio::geo::kilometers(80)});
    // Verify the number of results for each query is correct
    EXPECT_EQ(result[0].size(), 2);
    EXPECT_EQ(result[1].size(), 3);
    EXPECT_EQ(result[2].size(), 1);
}
