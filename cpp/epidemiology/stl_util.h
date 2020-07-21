#ifndef STL_UTIL_H
#define STL_UTIL_H

#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>
#include <cassert>

namespace epi
{

/**
 * @brief inserts element in a sorted vector, replacing items that are equal
 * precondition:  elements in the vector are partially sorted and unique according the predicate
 * postcondition: same as precondition, additionally contains exactly one element that is equal to item, 
 *                order of other items is conserved
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
    }
    else {
        vec.insert(lb, item);
    }
    return lb;
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
    using iterator       = typename IterPair::first_type;
    using const_iterator = typename IterPair::first_type;
    using value_type     = typename std::iterator_traits<iterator>::value_type;
    using reference      = typename std::iterator_traits<iterator>::reference;

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
 * template meta programming helper type 
 */
template <class... Ts>
using void_t = void;

/**
 * meta function to check type T for an existing stream output operator "<<"
 */
template <class T, class = void>
struct has_ostream_op : std::false_type {
};
template <class T>
struct has_ostream_op<T, void_t<decltype(std::declval<std::ostream&>() << std::declval<T>())>> : std::true_type {
};

/**
 * meta function to check type T for an existing equality comparison operator
 */
template <class T, class = void>
struct has_eq_op : std::false_type {
};
template <class T>
struct has_eq_op<T, void_t<decltype(std::declval<T>() == std::declval<T>())>> : std::true_type {
};

}

#endif //STL_UTIL_H