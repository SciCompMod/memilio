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

// empty struct to pass parameter packs
template <class... Ts>
struct typelist {
};

// function declaration used to remove Omittand from the type list of a tuple
template <class Omittand, class... Tags>
decltype(std::tuple_cat(std::declval<typename std::conditional<std::is_same<Omittand, Tags>::value, std::tuple<>,
                                                               std::tuple<Tags>>::type>()...))
    filter_tuple(std::tuple<Tags...>);

// function declaration used to replace type T by std::tuple
template <template <class...> class T, class... Args>
std::tuple<Args...> as_tuple(T<Args...>);

// function declaration used to replace std::tuple by type T
template <template <class...> class T, class... Args>
T<Args...> as_index(std::tuple<Args...>);

// remove all occurances of Omittand from the types in a std::tuple<types...>
template <class Omittand, class Tuple>
using filtered_tuple_t = decltype(filter_tuple<Omittand>(std::declval<Tuple>()));

// remove all occurances of Omittand from the types in an Index = IndexTemplate<types...>
template <class Omittand, template <class...> class IndexTemplate, class Index>
using filtered_index_t = decltype(as_index<IndexTemplate>(
    std::declval<filtered_tuple_t<Omittand, decltype(as_tuple(std::declval<Index>()))>>()));

} // namespace details

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