/* 
* Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
*
* Authors: Wadim Koslow, Daniel Abele
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
#ifndef PARAMETER_SPACE_H
#define PARAMETER_SPACE_H

#include "memilio/mobility/mobility.h"
#include "memilio/utils/memory.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/parameter_distributions.h"
#include "ode_secirvvs/model.h"

#include <assert.h>
#include <string>
#include <vector>
#include <random>
#include <memory>

namespace mio
{
namespace osecirvvs
{
    /**
     * draws a sample from the specified distributions for all parameters related to the demographics, e.g. population.
     * @param[inout] model Model including contact patterns for alle age groups
     */
    void draw_sample_demographics(Model& model);

    /**
     * draws a sample from the specified distributions for all parameters related to the infection.
     * @param[inout] model Model including contact patterns for alle age groups
     */
    void draw_sample_infection(Model& model);

    /** Draws a sample from Model parameter distributions and stores sample values
    * as Parameters parameter values (cf. UncertainValue and Parameters classes)
    * @param[inout] model Model including contact patterns for alle age groups
    */
    void draw_sample(Model& model);

    Graph<Model, MigrationParameters> draw_sample(Graph<Model, MigrationParameters>& graph, bool high);

} // namespace osecirvvs
} // namespace mio

#endif // PARAMETER_SPACE_H
