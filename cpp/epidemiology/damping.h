#include <vector>
#include <algorithm>

namespace epi
{

/**
 * This defined a damping factor for a
 * mitigation strategy for one point in time.
 */
class Damping
{
public:
    double day;
    double factor;

    Damping(double day_in, double factor_in);
};

template <typename T, typename Pred>
typename std::vector<T>::iterator insert_sorted(std::vector<T>& vec, T const& item, Pred pred)
{
    return vec.insert(std::upper_bound(vec.begin(), vec.end(), item, pred), item);
}

} // namespace epi