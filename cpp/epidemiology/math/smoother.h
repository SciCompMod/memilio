#include <cmath>

namespace epi
{

/**
 * @brief Returns the smoothed evaluation of a discrete jump of function values  
 * yleft and yright on xleft and xright, respectively, by using a cosine function.
 * 
 * @param x current evaluation point
 * @param xleft left boundary of independent variable
 * @param xright right boundary of independent variable
 * @param yleft function value at left boundary
 * @param yright function value at right boundary
 * @return double cosine-smoothed evaluation of discrete step function
 */
inline double smoother_cosine(double x, double xleft, double xright, double yleft, double yright)
{

    // f(x) = 0.5*(yleft - yright)*cos(pi/ descent_area*(x-day_upper_min))+0.5*(yleft + yright)
    return 0.5 * (yleft - yright) * std::cos(3.14159265358979323846 / (xright - xleft) * (x - xleft)) +
           0.5 * (yleft + yright);
}

} // namespace epi
