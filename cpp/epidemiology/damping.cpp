#include <epidemiology/damping.h>

#include <algorithm>
#include <cassert>

namespace
{

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

} // namespace

namespace epi
{

Damping::Damping(double day_in, double factor_in)
    : day(day_in)
    , factor(factor_in)
{
}

Dampings::Dampings()
    : m_dampings({{0.0, 1.0}})
{
}

void Dampings::add(const Damping& d)
{
    // make sure, the damping array is sorted
    insert_sorted_replace(m_dampings, d, [](const Damping& d1, const Damping& d2) {
        return d1.day < d2.day;
    });
}

double Dampings::get_factor(double day) const
{
    //returns the damping that is active at the specified time or the first damping
    assert(!m_dampings.empty());
    auto ub = std::upper_bound(begin(m_dampings), end(m_dampings), day, [](double d1, const Damping& d2) {
        return d1 < d2.day;
    });
    return (ub == begin(m_dampings) ? ub : ub - 1)->factor;
}

} // namespace epi
