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
    std::cout << "Forward AD" << std::endl;
    ad::gt1s<double>::type t1_x;
    ad::value(t1_x)             = 3.0;
    ad::derivative(t1_x)        = 1.0;
    ad::gt1s<double>::type t1_y = square(t1_x);
    std::cout << "value of square(" << ad::value(t1_x) << ") is " << ad::value(t1_y) << std::endl;
    std::cout << "forward derivative of square(" << ad::value(t1_x) << ") is " << ad::derivative(t1_y) << std::endl;

    std::cout << "Reverse AD" << std::endl;
    // create tape
    if (!ad::ga1s<double>::global_tape)
        ad::ga1s<double>::global_tape = ad::ga1s<double>::tape_t::create();
    ad::ga1s<double>::global_tape->reset();
    ad::ga1s<double>::type a1_x;
    ad::value(a1_x)      = 3.0;
    ad::derivative(a1_x) = 0.0;
    ad::ga1s<double>::global_tape->register_variable(a1_x);
    ad::ga1s<double>::type a1_y = square(a1_x);
    std::cout << "value of square(" << ad::value(a1_x) << ") is " << ad::value(t1_y) << std::endl;
    ad::ga1s<double>::global_tape->register_output_variable(a1_y);
    ad::derivative(a1_y) = 1.0;
    ad::ga1s<double>::global_tape->interpret_adjoint();
    std::cout << "adjoint derivative is " << ad::derivative(a1_x) << std::endl;
    ad::ga1s<double>::tape_t::remove(ad::ga1s<double>::global_tape);

    return 0;
}
