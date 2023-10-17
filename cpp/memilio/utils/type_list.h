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
#ifndef MIO_UTILS_TYPE_LIST_H_
#define MIO_UTILS_TYPE_LIST_H_

#include "memilio/io/io.h"
#include "memilio/utils/index.h"

#include <tuple>
#include <type_traits>
#include <utility>

namespace mio
{

template <size_t Index, class... Types>
struct type_at {
    template <template <class...> class TypeContainer>
    type_at(TypeContainer<Types...>)
    {
    }
    using type = typename std::tuple_element<Index, std::tuple<Types...>>::type;
    // std::remove_reference_t<decltype(std::get<Index>(std::declval<std::tuple<Types...>>()))>;
};

template <size_t Index, class... Types>
using type_at_t = typename type_at<Index, Types...>::type;

template <class Type, class... Types>
struct index_of {
    template <template <class...> class TypeContainer>
    index_of(TypeContainer<Types...>)
    {
    }

private:
    template <size_t Index = 0>
    static constexpr size_t index_of_impl()
    {
        if constexpr (std::is_same<Type, typename std::tuple_element<Index, std::tuple<Types...>>::type>::value) {
            return Index;
        }
        else {
            static_assert(Index < sizeof...(Types), "Type is not contained in given list.");
            return index_of_impl<Index + 1>();
        }
    }

public:
    static constexpr size_t value = index_of_impl();
};

template <class Type, class... Types>
constexpr size_t index_of_v = index_of<Type, Types...>::value;

/// @brief Collection of types. Each type is mapped to an index of type size_t.
template <class... Types>
struct TypeList {
    /**
     * @brief Get a type by index.
     * @tparam Index position of a type.
     * @return The type at position Index of TypeList.
     */
    template <size_t Index>
    constexpr auto get() const
    {
        return type_at_t<Index, Types...>();
    }

    /**
     * @brief Get index of a given type.
     * @tparam Type A type contained in TypeList.
     * @return Position of Type within TypeList.
     */
    template <class Type>
    constexpr size_t get() const
    {
        return index_of_v<Type, Types...>;
    }

    /// @brief returns the number of Types in TypeList
    constexpr size_t size() const
    {
        return sizeof...(Types);
    }
};

template <size_t Index, class... Types>
struct type_at<Index, TypeList<Types...>> : public type_at<Index, Types...> {
};

template <class Type, class... Types>
struct index_of<Type, TypeList<Types...>> : public index_of<Type, Types...> {
};

} // namespace mio

#endif // MIO_UTILS_TYPE_LIST_H_