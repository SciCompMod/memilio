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
#include "ide_sir_analytical/gregory_weights.h"

namespace mio
{
namespace isir
{
// Gregory weights corresponding to gregory_order where gregory_order corresponds to n0. The expected order of
// convergence is gregory_order+1.
std::vector<Eigen::MatrixX<ScalarType>> get_gregoryweights(size_t gregory_order)
{
    Eigen::MatrixX<ScalarType> gregoryWeights_sigma(gregory_order, gregory_order);
    Eigen::VectorX<ScalarType> gregoryWeights_omega(gregory_order + 1);
    ScalarType scale_gregory_weights;

    if (gregory_order == 1) {
        scale_gregory_weights = 2.;
        gregoryWeights_sigma << 1 / scale_gregory_weights;
        gregoryWeights_omega << 1. / scale_gregory_weights, 2. / scale_gregory_weights;
    }

    if (gregory_order == 2) {
        scale_gregory_weights = 12.;
        gregoryWeights_sigma << 5. / scale_gregory_weights, 14. / scale_gregory_weights, 5. / scale_gregory_weights,
            13. / scale_gregory_weights;
        gregoryWeights_omega << 5. / scale_gregory_weights, 13. / scale_gregory_weights, 12. / scale_gregory_weights;
    }

    // As in Pezzella.
    if (gregory_order == 3) {
        scale_gregory_weights = 24.;
        gregoryWeights_sigma << 9. / scale_gregory_weights, 27. / scale_gregory_weights, 27. / scale_gregory_weights,
            9. / scale_gregory_weights, 28. / scale_gregory_weights, 22. / scale_gregory_weights,
            9. / scale_gregory_weights, 28. / scale_gregory_weights, 23. / scale_gregory_weights;
        gregoryWeights_omega << 9. / scale_gregory_weights, 28. / scale_gregory_weights, 23. / scale_gregory_weights,
            24. / scale_gregory_weights;
    }

    // // Other weights.
    // if (gregory_order == 3) {
    //     scale_gregory_weights = 24.;
    //     gregoryWeights_sigma << 11. / scale_gregory_weights, 23. / scale_gregory_weights, 29. / scale_gregory_weights,
    //         11. / scale_gregory_weights, 24. / scale_gregory_weights, 24. / scale_gregory_weights,
    //         11. / scale_gregory_weights, 24. / scale_gregory_weights, 25. / scale_gregory_weights;
    //     gregoryWeights_omega << 9. / scale_gregory_weights, 28. / scale_gregory_weights, 23. / scale_gregory_weights,
    //         24. / scale_gregory_weights;
    // }

    // if (gregory_order == 4) {
    //     scale_gregory_weights = 720.;
    //     gregoryWeights_sigma << 311. / scale_gregory_weights, 796. / scale_gregory_weights,
    //         606. / scale_gregory_weights, 916. / scale_gregory_weights, //
    //         311. / scale_gregory_weights, 777. / scale_gregory_weights, 712. / scale_gregory_weights,
    //         652. / scale_gregory_weights, //
    //         311. / scale_gregory_weights, 777. / scale_gregory_weights, 693. / scale_gregory_weights,
    //         758. / scale_gregory_weights, //
    //         311. / scale_gregory_weights, 777. / scale_gregory_weights, 693. / scale_gregory_weights,
    //         739. / scale_gregory_weights, //
    //         gregoryWeights_omega << 251. / scale_gregory_weights, 897. / scale_gregory_weights,
    //         633. / scale_gregory_weights, 739. / scale_gregory_weights, 720. / scale_gregory_weights;
    // }

    return {gregoryWeights_sigma, gregoryWeights_omega};
}
} // namespace isir
} // namespace mio