/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Martin J. Kuehn
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
#include "memilio/io/mobility_io.h"

// wrapper function to print out matrix entries by gdb's 'print get_element(M,1,1)'
// (GDB doesn't support calling the overloaded operator())
double get_element(Eigen::MatrixXd const& m, int i, int j)
{
    return m(i, j);
}

int main()
{
    // Place text file needs to be in working directory build/examples/ and
    // start from within examples folder in build directory
    auto twitter_migration_2018 = mio::read_mobility_formatted("2018_lk_matrix.txt");
}
