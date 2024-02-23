/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

#ifndef IDE_SECIR_INITIALIZE_FROM_ODE_H
#define IDE_SECIR_INITIALIZE_FROM_ODE_H

#include "ide_secir/model_ide.h"
#include "ide_secir/parameters.h"
#include "ode_secir/model.h"
#include "infection_state.h"
#include "memilio/config.h"
#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"

namespace mio
{
namespace isecir
{

void compute_initial_flows_from_ode_compartments(mio::osecir::Model& model_ode, mio::isecir::Model& model_ide,
                                                 mio::TimeSeries<ScalarType> secihurd_ode, ScalarType t0_ide,
                                                 ScalarType dt);

} // namespace isecir
} // namespace mio

#endif //IDE_SECIR_INITIALIZE_FROM_ODE_H