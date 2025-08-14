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

#include "memilio/geography/tree.h"
#include "memilio/geography/locations.h"
#include "random_number_test.h"

using TestGeography = RandomNumberTest;

/**
 * @brief Test comparing geographical locations for equality.
 */
TEST_F(TestGeography, compareGeographicalLocation)
{
    // Set a geographical location for the location.
    mio::geo::GeographicalLocation geographical_location  = {10.5100470359749, 52.2672785559812};
    mio::geo::GeographicalLocation geographical_location2 = {10.5100470359749, 52.2672785559812};
    // Verify that the set geographical location matches the expected values.
    EXPECT_TRUE(geographical_location == geographical_location2);
}

/**
 * @brief Test comparing geographical locations for inequality.
 */
TEST_F(TestGeography, compareGeographicalLocation2)
{
    // Set a geographical location for the location.
    mio::geo::GeographicalLocation geographical_location  = {10.5100470359749, 52.2672785559812};
    mio::geo::GeographicalLocation geographical_location2 = {10.5100470309749, 52.2672785559812};
    // Verify that the set geographical location matches the expected values.
    EXPECT_FALSE(geographical_location == geographical_location2);
}

/**
 * @brief Test calculating the distance between two locations
 */
TEST_F(TestGeography, Distance)
{
    auto bonn   = mio::geo::GeographicalLocation(50.7333, 7.1000);
    auto berlin = mio::geo::GeographicalLocation(52.5200, 13.4050);
    EXPECT_DOUBLE_EQ(bonn.distance(berlin), berlin.distance(bonn));
    auto distance = 478.2;
    EXPECT_LT(abs(bonn.distance(berlin) - distance), 0.1);
}

/**
 * @brief Test the default r-Tree Constructor
 */
TEST_F(TestGeography, rtreeConstructionNoData)
{
    EXPECT_NO_THROW(mio::geo::RTree());
}

/**
 * @brief Test the r-Tree Constructor given a vector
 */
TEST_F(TestGeography, rtreeConstructionVector)
{
    std::vector<mio::geo::GeographicalLocation> locations;
    locations.push_back(mio::geo::GeographicalLocation(50.783, 6.083));
    locations.push_back(mio::geo::GeographicalLocation(52.083, 7.017));
    locations.push_back(mio::geo::GeographicalLocation(53.667, 10.233));
    EXPECT_NO_THROW(mio::geo::RTree tree(locations));
}

/**
 * @brief Test the r-Tree Constructor given a range
 */
TEST_F(TestGeography, rtreeConstructionRange)
{
    std::vector<mio::geo::GeographicalLocation> locations;
    locations.push_back(mio::geo::GeographicalLocation(50.783, 6.083));
    locations.push_back(mio::geo::GeographicalLocation(52.083, 7.017));
    locations.push_back(mio::geo::GeographicalLocation(53.667, 10.233));
    EXPECT_NO_THROW(mio::geo::RTree(locations.begin(), locations.end()));
}

/**
 * @brief Test the size function of an r-Tree
 */
TEST_F(TestGeography, rtreesize)
{
    std::vector<mio::geo::GeographicalLocation> locations;
    locations.push_back(mio::geo::GeographicalLocation(50.783, 6.083));
    locations.push_back(mio::geo::GeographicalLocation(52.083, 7.017));
    locations.push_back(mio::geo::GeographicalLocation(53.667, 10.233));
    auto rtree = mio::geo::RTree(locations.begin(), locations.end());
    EXPECT_EQ(rtree.size(), 3);
}

/**
 * @brief Test the nearest neighbours query of an r-Tree
 */
TEST_F(TestGeography, rtreeNN)
{
    std::vector<mio::geo::GeographicalLocation> locations;
    locations.push_back(mio::geo::GeographicalLocation(50.783, 6.083));
    locations.push_back(mio::geo::GeographicalLocation(52.083, 7.017));
    locations.push_back(mio::geo::GeographicalLocation(53.667, 10.233));
    auto rtree = mio::geo::RTree(locations.begin(), locations.end());
    EXPECT_EQ(rtree.nearest_neighbor_indices(mio::geo::GeographicalLocation(50.781, 6.080), 1)[0], 0);
}

/**
 * @brief Test the in-range query of an r-Tree
 */
TEST_F(TestGeography, rtreeinrange_approx)
{
    std::vector<mio::geo::GeographicalLocation> locations;
    locations.push_back(mio::geo::GeographicalLocation(50.783, 6.083));
    locations.push_back(mio::geo::GeographicalLocation(52.083, 7.017));
    locations.push_back(mio::geo::GeographicalLocation(53.667, 10.233));
    auto rtree = mio::geo::RTree(locations.begin(), locations.end());
    EXPECT_EQ(rtree.inrange_indices_approximate(mio::geo::GeographicalLocation(50.933501, 6.875124), 150).size(), 2);
}

/**
 * @brief Test the exact in-range query of an r-Tree
 */
TEST_F(TestGeography, rtreeinrange)
{
    std::vector<mio::geo::GeographicalLocation> locations;
    locations.push_back(mio::geo::GeographicalLocation(50.783, 6.083));
    locations.push_back(mio::geo::GeographicalLocation(52.083, 7.017));
    locations.push_back(mio::geo::GeographicalLocation(53.667, 10.233));
    auto rtree = mio::geo::RTree(locations.begin(), locations.end());
    EXPECT_EQ(rtree.inrange_indices(mio::geo::GeographicalLocation(50.933501, 6.875124), 150).size(), 2);
}

/**
 * @brief Test the in-range query of an r-Tree with multiple radii
 */
TEST_F(TestGeography, rtreeinrange_multiple_radii)
{
    std::vector<mio::geo::GeographicalLocation> locations;
    locations.push_back(mio::geo::GeographicalLocation(50.783, 6.083));
    locations.push_back(mio::geo::GeographicalLocation(52.083, 7.017));
    locations.push_back(mio::geo::GeographicalLocation(53.667, 10.233));
    auto rtree  = mio::geo::RTree(locations.begin(), locations.end());
    auto result = rtree.inrange_indices_query(mio::geo::GeographicalLocation(51.492599, 7.451810), {130, 310, 80});
    EXPECT_EQ(result[0].size(), 2);
    EXPECT_EQ(result[1].size(), 3);
    EXPECT_EQ(result[2].size(), 1);
}
