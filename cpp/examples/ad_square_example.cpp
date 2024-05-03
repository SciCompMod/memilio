/*
* Copyright (C) 2020-2024 MEmilio
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

template <typename FP>
FP square(FP x)
{
    return x * x;
}

int main()
{
    /**** Forward AD ****/
    std::cout << "Forward AD" << std::endl;
    ad::gt1s<double>::type t1_x; // tangent-linear AD type
    ad::value(t1_x)             = 3.0; // set value (this corresponds to the non-AD part)
    ad::derivative(t1_x)        = 1.0; // set tangent-linear derivative seed
    ad::gt1s<double>::type t1_y = square(t1_x);
    std::cout << "value of square(" << ad::value(t1_x) << ") is " << ad::value(t1_y) << std::endl;
    std::cout << "forward derivative of square(" << ad::value(t1_x) << ") is " << ad::derivative(t1_y) << std::endl;

    /**** Reverse AD ****/
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
