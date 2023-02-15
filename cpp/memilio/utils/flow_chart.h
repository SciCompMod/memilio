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

#include <tuple>
#include <type_traits>

namespace mio
{

namespace details
{

// base type to filter tuple Tags
template <class Omittand, class Tag, class... Tags>
struct FilteredTuple {
};
// 1st specialization: omitts Omittand, used by the 3rd generalization
template <class Omittand>
struct FilteredTuple<Omittand, std::tuple<Omittand>> {
    using Type = std::tuple<>;
};
// 2nd specialization: accept Tag, confirmed to not be of type Omittand
// (the 1st specialization is preferred, since it uses an explicit Type)
template <class Omittand, class Tag>
struct FilteredTuple<Omittand, std::tuple<Tag>> {
    using Type = std::tuple<Tag>;
};
// 3rd/main specialization: recursively uses tuple_cat to concatonate Tag(s) != Omittand
template <class Omittand, class Tag, class... Tags>
struct FilteredTuple<Omittand, std::tuple<Tag, Tags...>> {
    using Type = decltype(std::tuple_cat(std::declval<typename FilteredTuple<Omittand, std::tuple<Tag>>::Type>(),
                                         std::declval<typename FilteredTuple<Omittand, std::tuple<Tags...>>::Type>()));
};

template <class Omittand, class... Tags>
using filtered_tuple_t = typename FilteredTuple<Omittand, Tags...>::Type;

template <template <class...> class T, class... Args>
constexpr std::tuple<Args...> as_tuple(T<Args...>);

template <template <class...> class T, class... Args>
constexpr T<Args...> as_index(std::tuple<Args...>);

template <class Omittand, template <class...> class IndexTemplate, class Index>
using filtered_index_t = decltype(as_index<IndexTemplate>(
    std::declval<typename FilteredTuple<Omittand, decltype(as_tuple(std::declval<Index>()))>::Type>()));

} // namespace details

template <class Status, Status Source, Status Target>
struct Flow {
    using Type                 = Status;
    const static Status source = Source;
    const static Status target = Target;
};

template <class... Flows>
class FlowChart
{
public:
    template <size_t index>
    constexpr auto get() const
    {
        return std::get<index>(m_flows);
    }

    template <class Flow>
    constexpr size_t get() const
    {
        return get_impl<0>(Flow());
    }

    constexpr size_t size() const
    {
        return sizeof...(Flows);
    }

private:
    template <class... Types>
    struct typelist {
    };

    template <size_t index, class Flow>
    inline constexpr std::enable_if_t<!std::is_same<Flow, decltype(FlowChart<Flows...>().get<index>())>::value, size_t>
    get_impl(Flow f) const
    {
        static_assert(index < sizeof...(Flows), "In get_impl<Flow>(): Flow is not in Flows...");
        return get_impl<index + 1>(f);
    }

    template <size_t index, class Flow>
    inline constexpr std::enable_if_t<std::is_same<Flow, decltype(FlowChart<Flows...>().get<index>())>::value, size_t>
    get_impl(Flow) const
    {
        return index;
    }

    std::tuple<Flows...> m_flows;
};

} // namespace mio

#endif