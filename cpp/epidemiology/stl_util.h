#ifndef STL_UTIL_H
#define STL_UTIL_H

#include <vector>
#include <algorithm>
#include <utility>
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

}

#endif //STL_UTIL_H