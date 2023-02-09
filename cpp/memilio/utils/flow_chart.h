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

template <class Status, Status Source, Status Target>
struct Flow {
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
        return get_impl(typelist<Flow>(), typelist<>(), typelist<Flows...>());
    }

    constexpr size_t size() const
    {
        return sizeof...(Flows);
    }

private:
    template <class... Types>
    struct typelist {
    };

    template <class Flow, class... Indices>
    constexpr size_t get_impl(typelist<Flow> f, typelist<Indices...>, typelist<>) const
    {
        static_assert(sizeof...(Indices) != sizeof...(Flows), "ERROR IN get<Flow>() - could not find flow");
        static_assert(sizeof...(Indices) == sizeof...(Flows), "ERROR IN get<Flow>() - could not find flow");
        return get_impl(f, typelist<Indices..., int>());
    }

    template <class Flow, class... Indices, class F, class... Fs>
    inline constexpr std::enable_if_t<!std::is_same<Flow, F>::value, size_t>
    get_impl(typelist<Flow> f, typelist<Indices...>, typelist<F, Fs...>) const
    {
        return get_impl(f, typelist<Indices..., int>(), typelist<Fs...>());
    }

    template <class Flow, class... Indices, class F, class... Fs>
    inline constexpr std::enable_if_t<std::is_same<Flow, F>::value, size_t>
    get_impl(typelist<Flow>, typelist<Indices...>, typelist<F, Fs...>) const
    {
        return sizeof...(Indices);
    }

    std::tuple<Flows...> m_flows;
};

#endif