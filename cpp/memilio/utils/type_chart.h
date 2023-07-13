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
#ifndef MIO_UTILS_TYPE_CHART_H_
#define MIO_UTILS_TYPE_CHART_H_

#include "memilio/io/io.h"
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
     * @brief Get a type by index.
     * @tparam Index position of a type.
     * @return The type at position Index of TypeChart.
     */
    template <size_t Index>
    constexpr auto get() const
    {
        return std::get<Index>(m_types);
    }

    /**
     * @brief Get index of a given type.
     * @tparam Type A type contained in TypeChart.
     * @return Position of Type within TypeChart.
     */
    template <class Type>
    constexpr size_t get() const
    {
        return get_impl<0>(mio::Tag<Type>());
    }

    /// @brief returns the number of Types in TypeChart
    constexpr size_t size() const
    {
        return sizeof...(Types);
    }

private:
    /// @brief Iterates Index via recursion, as long as Type is not the Index-th position of Types.
    template <size_t Index, class Type>
    inline constexpr std::enable_if_t<
        !std::is_same<Type, typename std::tuple_element<Index, std::tuple<Types...>>::type>::value, size_t>
    get_impl(mio::Tag<Type> f) const
    {
        static_assert(Index < sizeof...(Types), "Type is not contained in TypeChart");
        return get_impl<Index + 1>(f);
    }

    /// @brief Returns Index. This ends the recursion of get_impl, when Type is the Index-th position of Types.
    template <size_t Index, class Type>
    inline constexpr std::enable_if_t<
        std::is_same<Type, typename std::tuple_element<Index, std::tuple<Types...>>::type>::value, size_t>
    get_impl(mio::Tag<Type>) const
    {
        return Index;
    }

    std::tuple<Types...> m_types; ///< Store types. Never instanced.
};

} // namespace mio

#endif // MIO_UTILS_TYPE_CHART_H_