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
#ifndef STL_UTIL_H
#define STL_UTIL_H

#include "memilio/utils/metaprogramming.h"

#include <vector>
#include <algorithm>
#include <utility>
#include <iterator>
#include <iostream>
#include <sstream>
#include <cstring>
#include <cassert>
#include <memory>

namespace mio
{

/**
 * @brief inserts element in a sorted vector, replacing items that are equal
 * precondition:  elements in the vector are partially sorted and unique according the predicate
 * postcondition: same as precondition, additionally contains exactly one element that is equal to item, 
 *                order of other items is preserved
 * @param vec vector where item will be inserted
 * @param item item to insert
 * @param pred binary comparator, pred(item, a) returns true if item should go before element a,
 *                                pred(a, item) returns true if element a should go before item
 * @return iterator to inserted or replaced item in vec
 */
template <typename T, typename Pred>
typename std::vector<T>::iterator insert_sorted_replace(std::vector<T>& vec, T const& item, Pred pred)
{
    auto bounds = std::equal_range(begin(vec), end(vec), item, pred);
    auto lb     = bounds.first;
    auto ub     = bounds.second;
    assert(ub - lb <= 1); //input vector contains at most one item that is equal to the new item
    if (ub - lb == 1) {
        *lb = item;
        return lb;
    }
    else {
        return vec.insert(lb, item);
    }
}

template <typename T>
typename std::vector<T>::iterator insert_sorted_replace(std::vector<T>& vec, T const& item)
{
    return insert_sorted_replace(vec, item, std::less<T>());
}

/**
 * @brief immutable random access range, e.g. a piece of a vector, represented by a pair of iterators
 * immutable means no elements can be added or removed. the elements themselves 
 * can be modified if the iterators that make up the range allow it.
 */
template <class IterPair>
class Range
{
public:
    using Iterators              = IterPair;
    using iterator               = typename IterPair::first_type;
    using const_iterator         = typename IterPair::first_type;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    using value_type             = typename std::iterator_traits<iterator>::value_type;
    using reference              = typename std::iterator_traits<iterator>::reference;

    Range(IterPair iter_pair)
        : m_iter_pair(iter_pair)
    {
    }

    /** @brief index operator.
     * constant complexity if random access iterator, linear otherwise
     */
    reference operator[](size_t idx) const
    {
        auto it = begin();
        std::advance(it, idx);
        return *it;
    }

    reference back() const
    {
        return *(--end());
    }

    size_t size() const
    {
        return static_cast<size_t>(std::distance(begin(), end()));
    }

    auto begin() const
    {
        return m_iter_pair.first;
    }

    auto end() const
    {
        return m_iter_pair.second;
    }

    template <class T = typename std::iterator_traits<iterator>::iterator_category,
              class   = std::enable_if_t<std::is_base_of<std::bidirectional_iterator_tag, T>::value>>
    auto rbegin() const
    {
        return std::make_reverse_iterator(end());
    }

    template <class T = typename std::iterator_traits<iterator>::iterator_category,
              class   = std::enable_if_t<std::is_base_of<std::bidirectional_iterator_tag, T>::value>>
    auto rend() const
    {
        return std::make_reverse_iterator(begin());
    }

private:
    IterPair m_iter_pair;
};

/**
 * @brief factories for template argument deduction
 */
template <class IterPair>
auto make_range(IterPair&& p)
{
    return Range<std::remove_reference_t<std::remove_cv_t<IterPair>>>{p};
}

template <class Iter1, class Iter2>
auto make_range(Iter1&& iter1, Iter2&& iter2)
{
    return make_range(std::make_pair(iter1, iter2));
}

/**
 * meta function to check type T for an existing stream output operator "<<"
 */
template <class T>
using ostream_op_expr_t = decltype(std::declval<std::ostream&>() << std::declval<T>());
template <class T>
using has_ostream_op = is_expression_valid<ostream_op_expr_t, T>;

/**
 * meta function to check type T for an existing equality comparison operator
 */
template <class T>
using eq_op_t = decltype(std::declval<T>() == std::declval<T>());
template <class T>
using has_eq_op = is_expression_valid<eq_op_t, T>;

/**
 * meta function to check type T for beeing an iterator
 */
template <class T>
using is_iterator_expr_t = typename std::iterator_traits<T>::iterator_category;
template <class T>
using is_iterator = is_expression_valid<is_iterator_expr_t, T>;

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

} // namespace mio

#endif //STL_UTIL_H
