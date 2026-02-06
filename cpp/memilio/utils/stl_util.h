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
#ifndef MIO_UTIL_STL_UTIL_H
#define MIO_UTIL_STL_UTIL_H

#include <array>
#include <numeric>
#include <ranges>
#include <vector>
#include <algorithm>
#include <utility>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include <cassert>
#include <memory>

namespace mio
{

/**
 * @brief Adds manipulators for width, (fixed) precision and fill character to an ostream.
 * Note that the formatting is consumed by the first output given to the ostream.
 * @param out Any std::ostream.
 * @param width Minimum width of the output.
 * @param precision The exact number of decimals (used only for numbers).
 * @param fill [Default: A space ' '] The character used for padding.
 * @return Returns a reference to out.
 */
inline std::ostream& set_ostream_format(std::ostream& out, size_t width, size_t precision, char fill = ' ')
{
    // Note: ostream& operator<< returns a reference to itself
    return out << std::setw(width) << std::fixed << std::setprecision(precision) << std::setfill(fill);
}

/**
 * @brief inserts element in a sorted vector, replacing items that are equal
 * precondition:  elements in the vector are partially sorted and unique according the predicate
 * postcondition: same as precondition, additionally contains exactly one element that is equal to item, 
 *                order of other items is preserved
 * @param vec vector where item will be inserted
 * @param item item to insert
 * @param pred binary comparator, pred(item, a) returns true if item should go before element a,
 *                                pred(a, item) returns true if element a should go before item
 */
template <typename T, typename Pred>
void insert_sorted_replace(std::vector<T>& vec, T const& item, Pred pred)
{
    auto bounds = std::equal_range(begin(vec), end(vec), item, pred);
    auto lb     = bounds.first;
    auto ub     = bounds.second;
    assert(ub - lb <= 1); //input vector contains at most one item that is equal to the new item
    if (ub - lb == 1) {
        *lb = item;
        return;
    }
    else {
        vec.insert(lb, item);
        return;
    }
}

template <typename T>
void insert_sorted_replace(std::vector<T>& vec, T const& item)
{
    return insert_sorted_replace(vec, item, std::less<T>());
}

/**
 * @brief A range of items defined by two iterators.
 * Models an immutable random access range, e.g. a piece of a vector. Immutable means that elements can not be added or
 * removed, while the elements themselves can be modified, if the iterators that make up the range allow it.
 * @tparam Iterator An iterator to the beginning of the range.
 * @tparam Sentinel An iterator to the end of the range. Defaults to the same type as Iterator.
 */
template <class Iterator, class Sentinel = Iterator>
class Range : public std::ranges::subrange<Iterator, Sentinel>
{
public:
    using Base = std::ranges::subrange<Iterator>;

    using iterator               = Iterator;
    using const_iterator         = iterator;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    using value_type             = typename std::iterator_traits<iterator>::value_type;
    using reference              = typename std::iterator_traits<iterator>::reference;

    /// @brief Directly construct a Range from two iterators.
    Range(Iterator begin, Sentinel end)
        : Base(std::move(begin), std::move(end))
    {
    }

    /// @brief Construct a Range from another range.
    Range(auto range)
        : Base(std::move(range).begin(), std::move(range).end())
    {
    }

    /// @brief Construct a Range from an std::pair of iterators.
    template <class I, class S>
    Range(std::pair<I, S> iterator_pair)
        : Base(std::move(iterator_pair).first, std::move(iterator_pair).second)
    {
    }

    decltype(auto) rbegin() const
        requires(std::is_base_of_v<std::bidirectional_iterator_tag,
                                   typename std::iterator_traits<iterator>::iterator_category>)
    {
        return std::make_reverse_iterator(Base::end());
    }

    decltype(auto) rend() const
        requires(std::is_base_of_v<std::bidirectional_iterator_tag,
                                   typename std::iterator_traits<iterator>::iterator_category>)
    {
        return std::make_reverse_iterator(Base::begin());
    }
};

/// @brief Deduction guide to correctly deduce range type when constructed from a pair.
template <class I, class S>
Range(std::pair<I, S> iterator_pair) -> Range<I, S>;

/// @brief Concept to check if type T has an existing stream output operator "<<".
template <class T>
concept HasOstreamOperator = requires(std::ostream os, T t) { os << t; };

namespace details
{
/**
 * length of a null terminated C string
 */
inline size_t string_length(const char* str)
{
    return std::strlen(str);
}

/**
 * length of a string (e.g std::string)
 */
template <class String>
size_t string_length(String&& str)
{
    return str.length();
}

/**
 * breaks the recursion of path_join_rec.
 */
inline void path_join_rec(std::stringstream&, bool)
{
}

/**
 * recursive template helper function to join paths
 * @param ss stream that collects the result
 * @param writeSeparator add separator before adding the next part of the path
 * @param head next part of the path to add
 * @param tail remaining parts of the path
 */
template <class Head, class... Tail>
void path_join_rec(std::stringstream& ss, bool writeSeparator, Head&& head, Tail&&... tail)
{
    if (writeSeparator) {
        ss << '/';
    }
    ss << head;
    path_join_rec(ss, string_length(head) > 0 && head[string_length(head) - 1] != '/', tail...);
}

} // namespace details

/** join one ore more strings with path separators.
 * Accepts mixed C strings or std::strings. 
 * 
 * example: 
 * 
 *     std::string hello("Hello");
 *     auto p = path_join(hello, "World"); //returns "Hello/World"
 * 
 * @param base first string
 * @param app zero or more other strings
 * @returns all inputs joined 
 */
template <class String, class... Strings>
std::string path_join(String&& base, Strings&&... app)
{
    std::stringstream ss;
    details::path_join_rec(ss, false, base, app...);
    auto path = ss.str();
    return path;
}

/**
 * converts a unique_ptr<T> to unique_ptr<U>.
 * behavior is similar to normal dynamic_cast except if the conversion 
 * is successful, the original unique_ptr<T> is now in a moved-from state and ownership
 * of the object has been transferred to the returned unique_ptr<U>.
 * @param base_ptr ptr to object to convert
 * @returns converted unique_ptr if object can be cast to U, default unique_ptr otherwise
 */
template <class U, class T>
std::unique_ptr<U> dynamic_unique_ptr_cast(std::unique_ptr<T>&& base_ptr)
{
    auto tmp_base_ptr = std::move(base_ptr);
    U* tmp            = dynamic_cast<U*>(tmp_base_ptr.get());
    std::unique_ptr<U> derived_ptr;
    if (tmp != nullptr) {
        tmp_base_ptr.release();
        derived_ptr.reset(tmp);
    }
    return derived_ptr;
}

/**
 * checks if there is an element in this range that matches a predicate
 */
template <class Iter, class Pred>
bool contains(Iter b, Iter e, Pred p)
{
    return find_if(b, e, p) != e;
}

/**
 * Get an std::array that contains all members of an enum class.
 * The enum class must be a valid index, i.e. members must be sequential starting at 0 and there must
 * be a member `Count` at the end, that will not be included in the array.
 * Example:
 * ```
 * enum class E { A, B, Count };
 * assert(enum_members<E>() == std::array<2, E>(E::A, E::B));
 * ``` 
 * @tparam T An enum class that is a valid index.
 * @return Array of all members of the enum class not including T::Count.
 */
template <class T>
constexpr std::array<T, size_t(T::Count)> enum_members()
{
    auto enum_members = std::array<T, size_t(T::Count)>{};
    auto indices      = std::array<std::underlying_type_t<T>, size_t(T::Count)>{};
    std::iota(indices.begin(), indices.end(), std::underlying_type_t<T>(0));
    std::transform(indices.begin(), indices.end(), enum_members.begin(), [](auto i) {
        return T(i);
    });
    return enum_members;
}

/**
 * @brief Defines generic Range type for IterPair of a vector.
 * When elements should be const, template argument accepts const type.
 * As vector is not defined for const values, the const specifier is removed and
 * the const_iterator is used. std::conditional tries to compile both cases, thus
 * we also need std::remove_const for the case where T is not const.
 */
template <class T>
using VectorRange =
    std::conditional_t<std::is_const_v<T>,
                       typename mio::Range<typename std::vector<std::remove_const_t<T>>::const_iterator>,
                       typename mio::Range<typename std::vector<std::remove_const_t<T>>::iterator>>;

} // namespace mio

#endif // MIO_UTIL_STL_UTIL_H
