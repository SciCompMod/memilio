#ifndef EPI_MATH_SMOOTHER_H
#define EPI_MATH_SMOOTHER_H

#include "epidemiology/utils/eigen.h"
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

/**
 * smoother_cosine as a matrix valued function.
 * @param x evaluation point
 * @param xleft left boundary x
 * @param xright right boundary x
 * @param yleft matrix expression, function value at left boundary
 * @param yright matrix expression, function value at right boundary
 * @return a matrix expression with yij = smoother_cosine(x, xleft, xright, yleftij, yrightij)
 */
template <class LeftExpr, class RightExpr,
          std::enable_if_t<std::is_base_of<Eigen::MatrixBase<std::decay_t<LeftExpr>>, std::decay_t<LeftExpr>>::value,
                           int> = 0>
auto smoother_cosine(double x, double xleft, double xright, const LeftExpr& yleft_expr, const RightExpr& yright_expr)
{
    return yleft_expr.binaryExpr(yright_expr, [=](auto yleft, auto yright) {
        return smoother_cosine(x, xleft, xright, yleft, yright);
    });
}

} // namespace epi

#endif //EPI_MATH_SMOOTHER_H