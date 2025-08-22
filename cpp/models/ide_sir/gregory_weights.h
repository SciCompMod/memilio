/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler
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

#ifndef IDESIR_GREGROY_WEIGHTS_H
#define IDESIR_GREGROY_WEIGHTS_H

#include "memilio/config.h"
#include "memilio/math/eigen.h"
// #include "memilio/math/eigen_util.h"
// #include "memilio/io/io.h"
// #include "memilio/utils/stl_util.h"
// #include "memilio/math/floating_point.h"

// #include <algorithm>
// #include <iterator>
// #include <vector>
// #include "memilio/utils/time_series.h"
// #include <Eigen/src/Core/util/Meta.h>
#include <cmath>

namespace mio
{
namespace isir
{
std::vector<Eigen::MatrixX<ScalarType>> get_gregoryweights(size_t gregory_order);
}
} // namespace mio
#endif // IDESIR_GREGROY_WEIGHTS_H