#include <cmath>

namespace epi
{

/**
 * @brief Returns the smoothed evaluation of a discrete jump of function values  
 * yleft and yright on xleft and xright, respectively, by using a cosine function.
 * If the input value is outside the given interval, yleft or yright are returned, respectively.
 *        { yleft,                                                                      for x <= xleft
 * f(x) = { yright,                                                                     for x >= xright
 *        { 0.5*(yleft - yright)*cos(pi/(xright-xleft)*(x-xleft))+0.5*(yleft + yright)  for x\in[xleft,xright]
 * @param x current evaluation point
 * @param xleft left boundary of independent variable
 * @param xright right boundary of independent variable
 * @param yleft function value at left boundary
 * @param yright function value at right boundary
 * @return double cosine-smoothed evaluation of discrete step function
 */
inline double smoother_cosine(double x, double xleft, double xright, double yleft, double yright)
{
    if (x <= xleft) {
        return yleft;
    }
    if (x >= xright) {
        return yright;
    }

    return 0.5 * (yleft - yright) * std::cos(3.14159265358979323846 / (xright - xleft) * (x - xleft)) +
           0.5 * (yleft + yright);
}

} // namespace epi