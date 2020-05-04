#include <epidemiology/damping.h>

namespace epi
{

Damping::Damping(double day_in, double factor_in)
    : day(day_in)
    , factor(factor_in)
{
}

} // namespace epi