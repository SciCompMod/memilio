/* 
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Martin J. Kuehn, Daniel Abele
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef MIO_MATH_SMOOTHER_H
#define MIO_MATH_SMOOTHER_H

#include "memilio/config.h"
#include "memilio/math/eigen.h"

#include <cmath>
#include <numbers>

namespace mio
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
 * @return Floating point cosine-smoothed evaluation of discrete step function
 */
template <typename FP>
inline FP smoother_cosine(FP x, FP xleft, FP xright, FP yleft, FP yright)
{
    using std::cos;

    if (x <= xleft) {
        return yleft;
    }
    if (x >= xright) {
        return yright;
    }

    return 0.5 * (yleft - yright) * cos(std::numbers::pi_v<ScalarType> / (xright - xleft) * (x - xleft)) +
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
template <typename FP, class LeftExpr, class RightExpr>
auto smoother_cosine(FP x, FP xleft, FP xright, const Eigen::MatrixBase<LeftExpr>& yleft_expr,
                     const Eigen::MatrixBase<RightExpr>& yright_expr)
{
    return yleft_expr.binaryExpr(yright_expr, [=](auto yleft, auto yright) {
        return smoother_cosine<FP>(x, xleft, xright, yleft, yright);
    });
}

} // namespace mio

#endif // MIO_MATH_SMOOTHER_H
