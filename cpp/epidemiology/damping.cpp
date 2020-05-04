#include <epidemiology/damping.h>

#include <algorithm>
#include <cassert>

namespace
{

template <typename T, typename Pred>
typename std::vector<T>::iterator insert_sorted(std::vector<T>& vec, T const& item, Pred pred)
{
    return vec.insert(std::upper_bound(vec.begin(), vec.end(), item, pred), item);
}

} // namespace

namespace epi
{

Damping::Damping(double day_in, double factor_in)
    : day(day_in)
    , factor(factor_in)
{
}

void Dampings::add(const Damping& d)
{
    // make sure, the damping array is sorted
    insert_sorted(m_dampings, d, [](const Damping& d1, const Damping& d2) {
        return d1.day < d2.day;
    });
}

double Dampings::get_factor(double day) const
{
    // we assume, that the data_array is ordered in ascending order
    size_t ilow  = 0;
    size_t ihigh = m_dampings.size() - 1;

    // check extrapolation cases
    if (day < m_dampings[ilow].day) {
        return m_dampings[ilow].factor;
    }

    if (day >= m_dampings[ihigh].day) {
        return m_dampings[ihigh].factor;
    }

    // now do the search
    while (ilow < ihigh - 1) {
        size_t imid = (ilow + ihigh) / 2;
        if (m_dampings[ilow].day <= day && day < m_dampings[imid].day) {
            ihigh = imid;
        }
        else if (m_dampings[imid].day <= day && day < m_dampings[ihigh].day) {
            ilow = imid;
        }
        else {
            // this case can only occur, if
            // input data are not ordered
            return 1e16;
        }
    }

    assert(m_dampings[ilow].day <= day && day < m_dampings[ilow + 1].day);
    return m_dampings[ilow].factor;
}

} // namespace epi
