/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Daniel Abele
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
#ifndef MIO_UTILS_METAPROGRAMMING_H
#define MIO_UTILS_METAPROGRAMMING_H

#include <type_traits>
#include <utility>

namespace ad
{
namespace internal
{
// Forward declaration of the AD type template
template <class Value, class Tape>
struct active_type;
} // namespace internal
} // namespace ad

namespace mio
{

namespace details
{
//non-copyable but trivally constructible and moveable type
struct NoCopy {
    NoCopy(const NoCopy&)            = delete;
    NoCopy()                         = default;
    NoCopy(NoCopy&&)                 = default;
    NoCopy& operator=(const NoCopy&) = delete;
    NoCopy& operator=(NoCopy&&)      = default;
};

//trivially constructible, copyable, moveable type
struct Empty {
};
} // namespace details

/// @brief Bool constant to check whether T is an integral type of 32 bits or less. Used for Json library.
template <class T>
using is_small_integral = std::integral_constant<bool, (std::is_integral_v<T> && sizeof(T) <= 4)>;
/// @brief Bool constant to check whether T is an integral type of 64 bits. Used for Json library.
template <class T>
using is_64bit_integral = std::integral_constant<bool, (std::is_integral_v<T> && sizeof(T) == 8)>;

/**
 * Defines a type that is not copy constructible or assignable if the specified condition is true.
 * Otherwise, defines a type that is copyable.
 * In both cases, the type will be trivially moveable and constructible.
 * To be used as a base type to make a class conditionally copyable.
 * Can be used to e.g. ensure that std::is_copy_constructible and std::is_copy_assignable is true only if 
 * the type can actually be copied (Note: this is not true for some STL containers, e.g., std::vector).
 */
template <bool Cond>
using not_copyable_if = std::conditional<Cond, details::NoCopy, details::Empty>;

/**
 * equivalent to not_copyable_if<Cond>::type. 
 * @see not_copyable_if
 */
template <bool Cond>
using not_copyable_if_t = typename not_copyable_if<Cond>::type;

namespace details
{

/**
 * @brief Get the type at position Index of list (Head, Tail...).
 * @tparam Index The index of the type to get.
 * @tparam Head, Tail The list to search in. May be empty.
 * @{
 */
template <std::size_t Index, class Head, class... Tail>
struct type_at_index_impl {
    using type = typename type_at_index_impl<Index - 1, Tail...>::type;
};
template <class Head, class... Tail>
struct type_at_index_impl<0, Head, Tail...> {
    using type = Head;
};
/** @} */

} // namespace details

/**
 * @brief Finds the type at the Index-th position in the list Types.
 * @tparam Index An index in `[0, sizeof...(Types))`.
 * @tparam Types A list of types.
 */
template <std::size_t Index, class... Types>
struct type_at_index {
    static_assert(Index < sizeof...(Types), "Index is too large for the list Types.");
    using type = typename details::type_at_index_impl<Index, Types...>::type;
};

/**
 * @brief The type at the Index-th position in the list Types.
 * Equivalent to type_at<Index, Types...>::type.
 * @see type_at_index.
 */
template <std::size_t Index, class... Types>
using type_at_index_t = typename type_at_index<Index, Types...>::type;

namespace details
{

/**
 * @brief Recursively searches the list (Head, Tail...) for Type.
 * @tparam Size Total size of the list to search.
 * @tparam Type to search.
 * @tparam Head, Tail The list to search in. May be empty.
 * @return The index of Type in the list (Head, Tail...), or the size of the list if Type is not in it.
 * @{
 */
template <std::size_t Size, class Type>
constexpr std::size_t index_of_impl()
{
    // this works both as an overload for empty lists as well as a "not found"
    return Size;
}
template <std::size_t Size, class Type, class Head, class... Tail>
constexpr std::size_t index_of_impl()
{
    // check if the type matches, otherwise call itself, omitting Head
    // this is significantly cheaper to compile compared to using an index and carrying the entire list,
    // as that needs additional work for looking up "Types[Index]", e.g. through type_at_index
    if constexpr (std::is_same_v<Type, Head>) {
        return Size - sizeof...(Tail) - 1;
    }
    else {
        return index_of_impl<Size, Type, Tail...>();
    }
}
/** @} */

} // namespace details

/**
 * Tests whether Type is in the list Types.
 * @tparam Type A type that may or may not be in Types.
 * @tparam Types A list of types.
 */
template <class Type, class... Types>
struct is_type_in_list
    : std::conditional_t<(details::index_of_impl<sizeof...(Types), Type, Types...>() < sizeof...(Types)),
                         std::true_type, std::false_type> {
};

/**
 * @brief Checks whether Type is in the list Types.
 * Equivalent to is_type_in_list<Type, Types...>::value.
 * @see is_type_in_list.
 */
template <class Type, class... Types>
constexpr bool is_type_in_list_v = is_type_in_list<Type, Types...>::value;

/**
 * @brief Finds the index of a Type in the list Types.
 * If Type does not have a unique index in Types, only the smallest is given as value.
 * @tparam Type A type contained in Types.
 * @tparam Types A list of types.
 */
template <class Type, class... Types>
struct index_of_type {
    static constexpr std::size_t value = details::index_of_impl<sizeof...(Types), Type, Types...>();
    static_assert(value < sizeof...(Types), "Type is not contained in given list.");
};

/**
 * @brief The index of Type in the list Types.
 * Equivalent to index_of_type<Type, Types...>::value.
 * @see index_of_type.
 */
template <class Type, class... Types>
constexpr std::size_t index_of_type_v = index_of_type<Type, Types...>::value;

/**
 * Tests whether the list Types contains any type multiple times.
 * @tparam Types A list of types.
 */
template <class... Types>
struct has_duplicates {
private:
    /**
     * @brief Checks if Types has a duplicate entry using an index sequence.
     * @tparam Indices Exactly the list '0, ... , sizeof...(Types) - 1'. Use std::make_index_sequence.
     * @return True if Types contains a duplicate type, false otherwise.
     */
    template <std::size_t... Indices>
    static constexpr bool has_duplicates_impl(std::index_sequence<Indices...>)
    {
        // index_of_type_v will always be equal to the index of the first occurance of a type,
        // while Indices contains its actual index. Hence, if there is any mismatch, then there is a duplicate.
        return ((index_of_type_v<Types, Types...> != Indices) || ...);
    }

public:
    static constexpr bool value = has_duplicates_impl(std::make_index_sequence<sizeof...(Types)>{});
};

/**
 * @brief Checks whether Type has any duplicates.
 * Equivalent to has_duplicates<Types...>::value.
 * @see is_type_in_list.
 */
template <class... Types>
constexpr bool has_duplicates_v = has_duplicates<Types...>::value;

/** 
 * @brief Detects whether a type is an automatic differentiation (AD) type. 
 * Intermediate result types from AD operations are not considered to be an AD type, as they use their own classes.
 * @{ 
 */
template <class T, class = void>
struct is_ad_type : public std::false_type {
};

template <class Value, class Tape>
struct is_ad_type<ad::internal::active_type<Value, Tape>> : public std::true_type {
};

template <class T>
constexpr bool is_ad_type_v = is_ad_type<T>::value;
/**@}*/

} // namespace mio

#endif // MIO_UTILS_METAPROGRAMMING_H
