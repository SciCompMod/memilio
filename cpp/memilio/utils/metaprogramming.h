/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#ifndef EPI_UTILS_METAPROGRAMMING_H
#define EPI_UTILS_METAPROGRAMMING_H

#include <type_traits>

namespace mio
{

namespace details
{
template <typename... Ts>
struct make_void {
    typedef void type;
};
} // namespace details

/**
 * utility for meta programming that produces void for any valid type.
 */
template <class... Ts>
using void_t = typename details::make_void<Ts...>::type;

namespace details
{
template <template <class...> class Expr, class X, class... T>
struct is_expression_valid : std::false_type {
};

template <template <class...> class Expr, class... T>
struct is_expression_valid<Expr, void_t<Expr<T...>>, T...> : std::true_type {
};
} // namespace details

/**
 * defines static constant value = true if Expr<T...> produces a valid type.
 * This technique is sometimes called detection idiom.
 * Can be used to detect e.g. a specific member functions:
 * @code
 *  struct Foo
 *  {
 *    void mem_fun(){}
 *  };
 * 
 *  template<class T>
 *  using mem_fun_t = decltype(std::declval<T>().mem_fun());
 * 
 *  static_assert(is_expression_valid<mem_fun_t, Foo>);
 * @endcode
 */
template <template <class...> class Expr, class... T>
struct is_expression_valid : details::is_expression_valid<Expr, void, T...> {
};

/**
 * Negates a type trait.
 * If B::value is true, then negation<B>::value is false.
 * @{
 */
template <class Trait>
struct negation : std::integral_constant<bool, (!Trait::value)> {
};
template <class Trait>
constexpr bool negation_v = negation<Trait>::value;
/**@}*/

/**
 * conjunction (logical and) of zero or more type traits with boolean values.
 * Does boolean shortcircuiting (like regular `if`) as expected. 
 * see https://en.cppreference.com/w/cpp/types/conjunction.
 * @{
 */
template <class...>
struct conjunction : std::true_type {
    //conjunction of no elements is true like in c++17
};
template <class B1>
struct conjunction<B1> : B1 {
    //conjunction of one element is identity
};
template <class B1, class... Bn>
struct conjunction<B1, Bn...> : std::conditional_t<bool(B1::value), conjunction<Bn...>, B1> {
    //conjunction of multiple elements is equal to the first element iff the first element is false.
    //otherwise its equal to the conjunction of the remaining elements.
};
template <class... Bs>
constexpr bool conjunction_v = conjunction<Bs...>::value;
/**@}*/

/**
* disjunction (logical or) of zero or more type traits with boolean values.
* Does boolean shortcircuiting (like regular `if`) as expected.
* see https://en.cppreference.com/w/cpp/types/disjunction.
* @{
*/
template<class...>
struct disjunction : std::false_type {
    //disjunction of no element is false like in c++17
};
template <class B1>
struct disjunction<B1> : B1 {
    //disjunction of one element is identity
};
template<class B1, class... Bn>
struct disjunction<B1, Bn...> : std::conditional<bool(B1::value), B1, disjunction<Bn...>> {
    //disjunction of mutliple elements is equal to the first element if the first element is true.
    //otherwise its equal to the disjunction of the remaining elements.
};
template<class... Bn>
constexpr bool disjunction_v = disjunction<Bn...>::value;
/**@}*/


namespace details
{
//non-copyable but trivally constructible and moveable type
struct NoCopy {
    NoCopy(const NoCopy&) = delete;
    NoCopy()              = default;
    NoCopy(NoCopy&&)      = default;
    NoCopy& operator=(const NoCopy&) = delete;
    NoCopy& operator=(NoCopy&&) = default;
};

//trivially constructible, copyable, moveable type
struct Empty {
};
} // namespace details

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

} // namespace mio

#endif
