/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: René Schmieding
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
#ifndef FLOW_CHART_H_
#define FLOW_CHART_H_

#include "memilio/utils/index.h"

#include <tuple>
#include <type_traits>
#include <utility>

namespace mio
{

namespace details
{

// helper struct to remove Omittand from the type list of a tuple
template <class Omittand, class Tuple>
struct FilteredTuple {
};
template <class Omittand>
struct FilteredTuple<Omittand, std::tuple<Omittand>> {
    // 1st specialization: omitts Omittand, used by the 3rd generalization
    using Type = std::tuple<>;
};
template <class Omittand, class Tag>
struct FilteredTuple<Omittand, std::tuple<Tag>> {
    // 2nd specialization: accept Tag, confirmed to not be of type Omittand
    // (the 1st specialization is preferred, since it uses an explicit Type)
    using Type = std::tuple<Tag>;
};
template <class Omittand, class Tag, class... Tags>
struct FilteredTuple<Omittand, std::tuple<Tag, Tags...>> {
    // 3rd/main specialization: recursively uses tuple_cat to concatonate Tag(s) != Omittand
    using Type = decltype(std::tuple_cat(std::declval<typename FilteredTuple<Omittand, std::tuple<Tag>>::Type>(),
                                         std::declval<typename FilteredTuple<Omittand, std::tuple<Tags...>>::Type>()));
};

// function declaration used to replace type T by std::tuple
template <template <class...> class T, class... Args>
std::tuple<Args...> as_tuple(T<Args...>);

// function declaration used to replace std::tuple by type T
template <template <class...> class T, class... Args>
T<Args...> as_index(std::tuple<Args...>);

// function declaration used to extract a Tag from type T
template <class Tag, template <class...> class T>
Tag get_tag(T<Tag>);

template <size_t I, class TagIndex, class Index, class Dim>
std::enable_if_t<(I == (TagIndex::size - 1) && TagIndex::size > 0), std::pair<size_t, size_t>>
flatten_index_by_tags(Index const& indices, Dim const& dimensions)
{
    using Tag = decltype(get_tag(mio::get<I>(std::declval<TagIndex>())));
    assert(get<Tag>(indices) < get<Tag>(dimensions));
    return {(size_t)mio::get<Tag>(indices), (size_t)mio::get<Tag>(dimensions)};
}

template <size_t I, class TagIndex, class Index, class Dim>
std::enable_if_t<(I<(TagIndex::size - 1) && TagIndex::size> 0), std::pair<size_t, size_t>>
flatten_index_by_tags(Index const& indices, Dim const& dimensions)
{
    using Tag = decltype(get_tag(mio::get<I>(std::declval<TagIndex>())));
    assert(mio::get<Tag>(indices) < mio::get<Tag>(dimensions));

    size_t val, prod;
    std::tie(val, prod) = flatten_index_by_tags<I + 1, TagIndex, Index, Dim>(indices, dimensions);

    return {val + (size_t)mio::get<Tag>(indices) * prod, prod * (size_t)mio::get<Tag>(dimensions)};
}

template <size_t I, class TagIndex, class Index, class Dim>
constexpr std::enable_if_t<(TagIndex::size == 0), std::pair<size_t, size_t>>
flatten_index_by_tags(Index const& /*indices*/, Dim const& /*dimensions*/)
{
    return {0, 0};
}

} // namespace details

// remove Omittand from the types in a std::tuple<types...>
template <class Omittand, class Tuple>
using filtered_tuple_t = typename details::FilteredTuple<Omittand, Tuple>::Type;

// remove Omittand from the types in an Index = IndexTemplate<types...>
template <class Omittand, template <class...> class IndexTemplate, class Index>
using filtered_index_t = decltype(details::as_index<IndexTemplate>(
    std::declval<filtered_tuple_t<Omittand, decltype(details::as_tuple(std::declval<Index>()))>>()));

/**
 * @brief calculates the flat index from a set of indices into a multidimensional array 
 *
 * @see flatten_index for more details
 *
 * @param indices a multiindex into a multidimensional array
 * @param dimensions a multiindex of sizes of each dimension of the array
 * @tparam TagIndex Index type, determining which Tags are used in which order for the flat index
 * @tparam Index Index type, containing (at least) all Tags in TagIndex
 * @tparam Dim Index type, containing (at least) all Tags in Index
 * @return the corresponding flat index
 */
template <class TagIndex, class Index, class Dim>
size_t flatten_index_by_tags(const Index& indices, const Dim& dimensions)
{
    return details::flatten_index_by_tags<0, TagIndex, Index, Dim>(indices, dimensions).first;
}

/// @brief A Flow defines a transition between two Compartments in a CompartmentalModel. Use with FlowChart
template <class Status, Status Source, Status Target>
struct Flow {
    using Type                 = Status;
    const static Status source = Source;
    const static Status target = Target;
};

/// @brief Collection of Flows defining the transitions between compartments in a CompartmentalModel
template <class... Flows>
class FlowChart
{
public:
    /**
     * @brief get a Flow by index
     * @tparam index position of a flow
     * @return Flow at position index of FlowChart
     */
    template <size_t index>
    constexpr auto get() const
    {
        return std::get<index>(m_flows);
    }

    /**
     * @brief get index of a given Flow
     * @tparam Flow a flow in FlowChart
     * @return position of Flow within FlowChart
     */
    template <class Flow>
    constexpr size_t get() const
    {
        return get_impl<0>(Flow());
    }

    /// @brief returns the number of Flows in FlowChart
    constexpr size_t size() const
    {
        return sizeof...(Flows);
    }

private:
    template <size_t index, class Flow>
    inline constexpr std::enable_if_t<
        !std::is_same<Flow, typename std::tuple_element<index, std::tuple<Flows...>>::type>::value, size_t>
    get_impl(Flow f) const
    {
        static_assert(index < sizeof...(Flows), "Flow is not contained in FlowChart");
        return get_impl<index + 1>(f);
    }

    template <size_t index, class Flow>
    inline constexpr std::enable_if_t<
        std::is_same<Flow, typename std::tuple_element<index, std::tuple<Flows...>>::type>::value, size_t>
    get_impl(Flow) const
    {
        return index;
    }

    std::tuple<Flows...> m_flows;
};

} // namespace mio

#endif