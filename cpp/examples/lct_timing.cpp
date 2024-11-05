/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Lena Ploetzke
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
#include "memilio/config.h"
#include "lct_timing.h"
#include <iostream>

int main(int argc, char** argv)
{
    const ScalarType tmax = 20;
    size_t warm_up        = 10;
    size_t num_runs       = 100;
    if (argc > 2) {
        warm_up  = std::stod(argv[1]);
        num_runs = std::stod(argv[2]);
    }
    constexpr size_t num_subcompartments = MY_COMPILE_TIME_VALUE;
    simulate_1<num_subcompartments>(warm_up, num_runs, tmax);
    return 0;
}
