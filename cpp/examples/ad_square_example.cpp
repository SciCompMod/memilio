/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Ralf Hannemann-Tamas
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

#include "ad/ad.hpp"
#include <iostream>

/*  This example computes the derivative of f(x) = x^2 in two different ways:
 *  1. Using the forward mode of algorithmic differentiation (AD)
 *  2. Using the reverse mode of algorithmic differentiation
 *
 *  Assume, a subroutine is given that computes y = f(x), where x is an n-dimenstional input vector
 *  and y is an m-dimensional output vector.
 *
 *  Forward mode AD (also known as tangent-linear mode)
 *  computes Jacobian-vector products of f(x). I.e., if J(x) denotes the Jacobian of f(x) and
 *  v is an n-dimensional seed vector, the forward mode computes the matrix-vector product
 *  J(x) * v. To this end, The data type of x is set to a forward mode AD data type that additionally to the
 *  original value, also stores derivative information. The original value is accessed by the command
 *  ad::value(x). The derivative information is accessed by ad::derivative(x) which is accordingly
 *  used to set the seed vector v, i.e., ad::derivative(x) = v.
 *  After the conmputation the Jacobian-vector product is accessed by the command ad::derivative(y).
 *
 *  Reverse mode (also known as adjoint mode, or back progpagation) computes vector-Jacobian products of f(x).
 *  I.e, if J(x) denotes the Jacobian of f(x) and
 *  w is an m-dimensional seed vector (which is a row vector), the reverse mode computes the vector-matrix product
 *  w * J(x). The data type of x is set to a reverse mode AD data type that additionally to the
 *  original value, also stores derivative information. The original value is accessed by the command
 *  ad::value(x). The derivative information is accessed by ad::derivative(x).
 *  While the forward mode computes derivatives in one sweep, the reverse mode computes is divided
 *  into two sweeps. A forward sweep for the computation, where intermediate results are stored on a temporary
 *  memory location called tape. After
 *  the forward sweep, the seed of the computation is set by the command ad::derivative(y) = w.
 *  And a reverse sweep, where the actual derivative computation takes place, has to be performed. The derivative information of
 *  x has to be initalized to zero (ad::derivative(x) = 0), before the tape can be interpreted backwards. Then,
 *  the tape can interpreted backwards
 *  and the vector-Jacobian product can be accessed the command ad::derivative(x).
 *
 *  For a concise introduction to AD we recommend the article:
 *  Verma, Arun. "An introduction to automatic differentiation." Current Science (2000): 804-807.
 */

/**
 * @brief square computes the square of the input parameter
 * @param x input parameter
 * @return the square of x.
 */
template <typename FP>
FP square(FP x)
{
    return x * x;
}

int main()
{
    /**** 1. Forward AD ****/
    std::cout << "Forward AD" << std::endl;
    ad::gt1s<double>::type t1_x; // tangent-linear AD type
    ad::value(t1_x)             = 3.0; // set value (this corresponds to the non-AD part)
    ad::derivative(t1_x)        = 1.0; // set tangent-linear derivative seed
    ad::gt1s<double>::type t1_y = square(t1_x);
    std::cout << "value of square(" << ad::value(t1_x) << ") is " << ad::value(t1_y) << std::endl;
    std::cout << "forward derivative of square(" << ad::value(t1_x) << ") is " << ad::derivative(t1_y) << std::endl;

    /**** 2. Reverse AD ****/
    std::cout << "Reverse AD" << std::endl;
    // create tape (allocation)
    if (!ad::ga1s<double>::global_tape)
        ad::ga1s<double>::global_tape = ad::ga1s<double>::tape_t::create();
    // clear tape
    ad::ga1s<double>::global_tape->reset();

    ad::ga1s<double>::type a1_x; // reverse-mode AD type
    ad::value(a1_x)      = 3.0; // set value (this corresponds to the non-AD part)
    ad::derivative(a1_x) = 0.0; // for reverse-mode the seed has to be set to zero
    ad::ga1s<double>::global_tape->register_variable(a1_x); // register input variable
    ad::ga1s<double>::type a1_y = square(a1_x);
    std::cout << "value of square(" << ad::value(a1_x) << ") is " << ad::value(a1_y) << std::endl;
    ad::ga1s<double>::global_tape->register_output_variable(a1_y);
    ad::derivative(a1_y) = 1.0;
    ad::ga1s<double>::global_tape
        ->interpret_adjoint(); // compute reverse-mode derivatives by evaluatin the tape backwards
    std::cout << "adjoint derivative is " << ad::derivative(a1_x) << std::endl; // access reverse-derivatives
    ad::ga1s<double>::tape_t::remove(ad::ga1s<double>::global_tape); // deallocate tape

    return 0;
}
