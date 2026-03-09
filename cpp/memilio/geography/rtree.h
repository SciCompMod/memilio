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
#ifndef MIO_GEOGRAPHY_RTREE_H
#define MIO_GEOGRAPHY_RTREE_H

#include "memilio/geography/distance.h"
#include "memilio/utils/back_inserter_second_element.h"

#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/strategies/cartesian/buffer_side_straight.hpp>
#include <boost/geometry/strategies/geographic/buffer_point_circle.hpp>

#include <cmath>
#include <iterator>
#include <type_traits>
#include <vector>

namespace mio
{
namespace geo
{
/**
 * @brief Concept for spherical location types.
 * 
 * The location type needs to provide get_latitude() and get_longitude() methods returning floating point types.
 */
template <class Location>
concept IsSphericalLocation = requires(const Location& loc) {
    std::is_floating_point_v<decltype(loc.get_latitude())> && std::is_floating_point_v<decltype(loc.get_longitude())>;
};

/**
 * @brief R-tree for spatial queries of geographical locations on the sphere.
 * 
 * Data structure to store spatial indices and allow for efficient in-range and nearest neighbour queries. 
 * Wraps the Boost::geometry::index::rtree. Can be initialized with a vector of geographical location data or a range.
 * The provided location data needs to provide get_latitude() and get_longitude().
 * The tree is initialised with a maximum number of elements per tree node of 16, which can be changed for different use 
 * cases. 
 */
class RTree
{
public:
    /**
     * @brief Construct a new RTree object without data.
     * 
     */
    RTree() = default;

    /**
     * @brief Construct a new RTree object with data given in a vector.
     * 
     * @param locations A vector of geographical locations, they need to provide get_latitude() and get_longitude().
     */
    template <IsSphericalLocation Location>
    RTree(const std::vector<Location>& locations)
        : m_rtree{}
    {
        for (size_t index = 0; index < locations.size(); index++) {
            Point point(locations[index].get_longitude(), locations[index].get_latitude());
            m_rtree.insert(Node(point, index));
        }
    }

    /**
     * @brief Return the number of data points stored in the RTree.
     * 
     * @return size of the tree.
     */
    auto size() const
    {
        return m_rtree.size();
    }

    /**
     * @brief Return the indices of the k nearest neighbors (i.e. nodes with the least distance) of a given location.
     * 
     * @param location Midpoint for the query, provides get_latitude() and get_longitude().
     * @param number The number of nearest neighbours to find.
     * @return Vector with indices of the nearest neighbours.
     */
    std::vector<size_t> nearest_neighbor_indices(const IsSphericalLocation auto& location, size_t number) const
    {
        using boost::geometry::index::nearest;
        using boost::geometry::index::query;

        Point point(location.get_longitude(), location.get_latitude());
        std::vector<size_t> indices;
        query(m_rtree, nearest(point, number), mio::back_inserter_second_element<std::vector<size_t>, Node>{indices});
        return indices;
    }

    /**
     * @brief Return the indices of the points within a given radius (in kilometers) where the circle is approximated by a polygon.
     * 
     * @param location Midpoint for the query, provides get_latitude() and get_longitude().
     * @param radius The radius of the query.
     * @return Vector with indices of the points found.
     */
    std::vector<size_t> in_range_indices_approximate(const IsSphericalLocation auto& location, Distance radius) const
    {
        using boost::geometry::index::covered_by;
        using boost::geometry::index::query;

        auto circle = create_circle_approximation(location, radius);
        std::vector<size_t> indices;
        query(m_rtree, covered_by(circle), mio::back_inserter_second_element<std::vector<size_t>, Node>{indices});
        return indices;
    }

    /**
     * @brief Return the indices of the points within each given radius for multiple radii at the same time.
     *
     * @param location Midpoint for the query, provides get_latitude() and get_longitude().
     * @param radii Vector containing the radii of the query.
     * @return Vector of vectors with indices of the points found.
     */
    std::vector<std::vector<size_t>> in_range_indices_query(const IsSphericalLocation auto& location,
                                                            const std::vector<Distance>& radii) const
    {
        using boost::geometry::distance;
        using boost::geometry::index::covered_by;
        using boost::geometry::index::query;

        auto max_radius = std::max_element(radii.begin(), radii.end());
        auto circle     = create_circle_approximation(location, *max_radius);
        std::vector<Node> nodes;
        query(m_rtree, covered_by(circle), std::back_inserter(nodes));
        auto midpoint = Point(location.get_longitude(), location.get_latitude());
        std::vector<std::vector<size_t>> result;
        for (const auto& radius : radii) {
            auto radius_in_meter = radius.meters();
            std::vector<size_t> indices;
            indices.reserve(nodes.size());
            for (auto& node : nodes) {
                if (distance(midpoint, node.first) < radius_in_meter) {
                    indices.push_back(node.second);
                }
            }
            result.push_back(indices);
        }
        return result;
    }

    /**
     * @brief Return the indices of the points within a given radius.
     * 
     * @param location Midpoint for the query, provides get_latitude() and get_longitude().
     * @param radius The radius of the query.
     * @return Vector with indices of the points found.
     *
     * Basically the same as \ref in_range_indices_approximate, but filters the result to make sure the points are within the radius.
     */
    std::vector<size_t> in_range_indices(const IsSphericalLocation auto& location, Distance radius) const
    {
        using boost::geometry::distance;
        using boost::geometry::index::covered_by;
        using boost::geometry::index::query;

        auto circle = create_circle_approximation(location, radius);

        std::vector<Node> nodes;
        query(m_rtree, covered_by(circle), std::back_inserter(nodes));
        auto midpoint = Point(location.get_longitude(), location.get_latitude());
        std::vector<size_t> indices;
        indices.reserve(nodes.size());
        for (auto& node : nodes) {
            if (distance(midpoint, node.first) < radius.meters()) {
                indices.push_back(node.second);
            }
        }
        return indices;
    }

private:
    /// @brief Point stores longitude and latitude in this order.
    using Point = boost::geometry::model::point<double, 2, boost::geometry::cs::geographic<boost::geometry::degree>>;
    /// @brief Node stores a point and its associated index.
    using Node = std::pair<Point, size_t>;
    /// @brief MultiPolygon type for circle approximation.
    using MultiPolygon = boost::geometry::model::multi_polygon<boost::geometry::model::polygon<Point>>;
    /**
     * @brief Create a circle approximation object.
     *
     * @param location Midpoint, needs to provide get_latitude() and get_longitude().
     * @param radius
     * @param approximation_points Number of points used to approximate the circle. Default is 36, i.e. we build a 36-gon.
     * @return multi_polygon.
     */
    MultiPolygon create_circle_approximation(const IsSphericalLocation auto& location, Distance radius,
                                             size_t approximation_points = 36) const
    {
        using namespace boost::geometry::strategy::buffer;

        geographic_point_circle<> point_strategy(approximation_points);
        distance_symmetric<double> distance_strategy(radius.meters());
        join_round join_strategy;
        end_round end_strategy;
        side_straight side_strategy;

        Point midpoint(location.get_longitude(), location.get_latitude());

        MultiPolygon circle;
        boost::geometry::buffer(midpoint, circle, distance_strategy, side_strategy, join_strategy, end_strategy,
                                point_strategy);
        return circle;
    }

    // Arbitrarily chosen value - can be changed for better performance in different use cases.
    constexpr static size_t m_max_number_of_elements_per_tree_node = 16;
    boost::geometry::index::rtree<Node, boost::geometry::index::rstar<m_max_number_of_elements_per_tree_node>> m_rtree;
};

} // namespace geo
} // namespace mio

#endif // MIO_GEOGRAPHY_RTREE_H
