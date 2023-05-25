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
#ifndef TYPE_CHART_H_
#define TYPE_CHART_H_

#include "memilio/utils/index.h"

#include <tuple>
#include <type_traits>
#include <utility>

namespace mio
{

/// @brief Collection of types. Each type is mapped to an index of type size_t.
template <class... Types>
class TypeChart
{
public:
    /**
     * @brief get a type by index
     * @tparam index position of a type
     * @return the type at position index of TypeChart
     */
    template <size_t index>
    constexpr auto get() const
    {
        return std::get<index>(m_types);
    }

    /**
     * @brief get index of a given type
     * @tparam Type a type in TypeChart
     * @return position of Type within TypeChart
     */
    template <class Type>
    constexpr size_t get() const
    {
        return get_impl<0>(Type());
    }

    /// @brief returns the number of Types in TypeChart
    constexpr size_t size() const
    {
        return sizeof...(Types);
    }

private:
    template <size_t index, class Type>
    inline constexpr std::enable_if_t<
        !std::is_same<Type, typename std::tuple_element<index, std::tuple<Types...>>::type>::value, size_t>
    get_impl(Type f) const
    {
        static_assert(index < sizeof...(Types), "Type is not contained in TypeChart");
        return get_impl<index + 1>(f);
    }

    template <size_t index, class Type>
    inline constexpr std::enable_if_t<
        std::is_same<Type, typename std::tuple_element<index, std::tuple<Types...>>::type>::value, size_t>
    get_impl(Type) const
    {
        return index;
    }

    std::tuple<Types...> m_types;
};

} // namespace mio

#endif // TYPE_CHART_H_