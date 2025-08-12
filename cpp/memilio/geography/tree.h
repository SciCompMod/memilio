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

#include "math.h"
#include "locations.h"
#include <boost/geometry/geometries/multi_polygon.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/strategies/cartesian/buffer_side_straight.hpp>
#include <boost/geometry/strategies/geographic/buffer_point_circle.hpp>
#include <concepts>

#include "boost/geometry/geometry.hpp"

#include <iterator>
#include <type_traits>
#include <vector>

namespace bg   = boost::geometry;
namespace bgi  = bg::index;
namespace bgsb = bg::strategy::buffer;

namespace mio
{
namespace geo
{
using Point = bg::model::point<double, 2, boost::geometry::cs::geographic<boost::geometry::degree>>;
typedef std::pair<Point, size_t> Node;

template <class Location>
concept IsSphericalLocation = requires(const Location& loc) {
    std::is_floating_point_v<decltype(loc.get_latitude())> && std::is_floating_point_v<decltype(loc.get_longitude())>;
};

template <class Iter>
concept IsSphericalLocationIterator = std::input_iterator<Iter> && IsSphericalLocation<decltype(*std::declval<Iter>())>;

/**
 * @brief R-Tree for spatial queries of geographical locations on the sphere
 * 
 * Wraps the Boost::geometry::index::rtree. Can be initialized with a vector of geographical location data or a range.
 * The provided location data needs to provide get_latitude() and get_longitude().
 */

class RTree
{
public:
    /**
     * @brief Construct a new RTree object without data
     * 
     */
    RTree() = default;

    /**
     * @brief Construct a new RTree object with data given in a vector
     * 
     * @param locations A vector of geographical locations, they need to provide get_latitude() and get_longitude().
     */
    template <IsSphericalLocation Location>
    RTree(const std::vector<Location>& locations)
        : rtree{}
    {
        for (size_t index = 0; index < locations.size(); index++) {
            Point point(locations[index].get_longitude(), locations[index].get_latitude());
            rtree.insert(Node(point, index));
        }
    }

    /**
     * @brief Construct a new RTree object with data given in a range
     * 
     * @param first The beginning of the range
     * @param last The end of the range
     * The provided location data needs to provide get_latitude() and get_longitude().
     */
    template <IsSphericalLocationIterator Iter>
    RTree(Iter first, Iter last)
        : rtree{}
    {
        size_t index = 0;
        while (first != last) {
            Point point(first->get_longitude(), first->get_latitude());
            rtree.insert(Node(point, index));
            ++first;
            ++index;
        }
    }
    /**
     * @brief Return the number of data points stored in the RTree
     * 
     * @return size of the tree 
     */
    auto size() const
    {
        return rtree.size();
    }

    /**
     * @brief Return the indices of the k nearest neighbours of a given location
     * 
     * @param location Midpoint for the query, provides get_latitude() and get_longitude()
     * @param number The number of nearest neighbours to find
     * @return Vector with indices of the nearest neighbours
     */
    auto nearest_neighbor_indices(const IsSphericalLocation auto& location, size_t number) const
    {
        Point point(location.get_longitude(), location.get_latitude());
        std::vector<size_t> indices;
        bgi::query(rtree, bgi::nearest(point, number), back_inserter_second_element<std::vector<size_t>>(indices));
        return indices;
    }

    /**
     * @brief Return the indices of the points within a given radius where the circle is approximated by a polygon
     * 
     * @param location Midpoint for the query, provides get_latitude() and get_longitude()
     * @param radius The radius of the query
     * @return Vector with indices of the points found
     */
    auto inrange_indices_approximate(const IsSphericalLocation auto& location, double radius) const
    {
        auto radius_in_meter = 1000 * radius;
        auto circle          = create_circle_approximation(location, radius_in_meter);
        std::vector<size_t> indices;
        bgi::query(rtree, bgi::covered_by(circle), back_inserter_second_element<std::vector<size_t>>(indices));
        return indices;
    }

    /**
     * @brief Return the indices of the points within a given radius
     * 
     * @param location Midpoint for the query, provides get_latitude() and get_longitude()
     * @param radius The radius of the query
     * @return Vector with indices of the points found
     *
     * Basically the same as \ref inrange_indices_approximate, but filters the result to make sure the points are within the radius.
     */
    auto inrange_indices(const IsSphericalLocation auto& location, double radius) const
    {
        auto radius_in_meter = 1000 * radius;

        auto circle = create_circle_approximation(location, radius_in_meter);
        std::vector<Node> nodes;
        bgi::query(rtree, bgi::covered_by(circle), std::back_inserter(nodes));
        auto midpoint = Point(location.get_longitude(), location.get_latitude());
        std::vector<size_t> indices;
        for (auto& node : nodes) {
            mio::log_info("Distance: {}", bg::distance(midpoint, node.first));
            if (bg::distance(midpoint, node.first) < radius_in_meter) {
                indices.push_back(node.second);
            }
        }
        return indices;
    }

private:
    /**
     * @brief Create a circle approximation object
     *
     * @param location Midpoint, needs to provide get_latitude() and get_longitude()
     * @param radius in meters
     * @return auto 
     */
    auto create_circle_approximation(const IsSphericalLocation auto& location, double radius) const
    {
        bgsb::geographic_point_circle<> point_strategy(36);
        bgsb::distance_symmetric<double> distance_strategy(radius);
        bgsb::join_round join_strategy;
        bgsb::end_round end_strategy;
        bgsb::side_straight side_strategy;

        Point midpoint(location.get_longitude(), location.get_latitude());

        bg::model::multi_polygon<bg::model::polygon<Point>> circle;
        bg::buffer(midpoint, circle, distance_strategy, side_strategy, join_strategy, end_strategy, point_strategy);
        return circle;
    }

    /**
     * @brief Back inserter that ignores the first element of pairs given to it
     */
    template <class Container>
    struct back_inserter_second_element {
        Container& container;
        back_inserter_second_element& operator*()
        {
            return *this;
        }
        back_inserter_second_element& operator++()
        {
            return *this;
        }
        back_inserter_second_element operator++(int)
        {
            return *this;
        }
        back_inserter_second_element operator=(const Node& node)
        {
            container.push_back(node.second);
            return *this;
        }
    };
    bgi::rtree<Node, bgi::rstar<16>> rtree;
};

} // namespace geo
} // namespace mio
