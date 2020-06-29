#include <epidemiology/damping.h>

#include <algorithm>
#include <epidemiology/stl_util.h>
#include <cassert>
#include <stdio.h>
#include <cmath>

namespace
{



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
    , m_smoothing(true)
{
}

std::vector<Damping> const& Dampings::get_dampings_vector() const
{
    return m_dampings;
}

void Dampings::set_smoothing(bool smoothing)
{
    m_smoothing = smoothing;
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

    // ub is the damping
    if (ub == begin(m_dampings)) { // only if day input is negative... (should not be the case...)
        return ub->factor;
    }
    else {
        if (m_smoothing) {
            double descent_area = 1.0; // smoothing of damping of 1 day max
            double day_upper    = 1e100; // large number; larger than maximum of days to be simulated
            if (ub < end(m_dampings)) { // if this is the case, there are at least two dampings in the list before

                day_upper    = (ub)->day;
                descent_area = std::min(1.0, day_upper - (ub - 1)->day);
            }

            if (day < day_upper - descent_area) {
                // here, the factor of the actual damping is return
                // printf("\n\n standard.. day %.6f factor %.6f ", day, (ub - 1)->factor);
                return (ub - 1)->factor;
            }
            else {
                // here, the transition of the actual to the next damping is smoothed
                // scale the cosine function such that (0,cos(0)) is mapped to (day_lower, fac_lower)
                // and (pi,cos(pi)) to (day_upper, fac_upper)

                double day_upper_min = day_upper - descent_area;
                double fac_lower     = (ub - 1)->factor; // factor at lower bound day
                double fac_upper     = (ub)->factor; // factor at upper bound day
                double ret;
                if (fac_lower >= fac_upper) {
                    // f(x) = 0.5*(fac_lower - fac_upper)*cos(pi/ descent_area*(x-day_upper_min))+0.5*(fac_lower + fac_upper)
                    ret = 0.5 * (fac_lower - fac_upper) *
                              std::cos(3.14159265358979323846 / descent_area * (day - day_upper_min)) +
                          0.5 * (fac_lower + fac_upper);
                }
                else {
                    ret = 0.5 * (fac_upper - fac_lower) *
                              std::cos(3.14159265358979323846 / descent_area * (descent_area + (day - day_upper_min))) +
                          0.5 * (fac_lower + fac_upper);
                }

                // printf("\n\n in descent.. day %.6f factor %.6f ", day, ret);

                return ret;
            }
        }
        else {
            return (ub - 1)->factor;
        }
    }
}

} // namespace epi
