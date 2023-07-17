/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kühn
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
#ifndef ODESECIRVVS_PARAMETER_SPACE_H
#define ODESECIRVVS_PARAMETER_SPACE_H

#include "memilio/mobility/metapopulation_mobility_instant.h"
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

/**
    * Draws samples for each model node in a graph.
    * Some parameters are shared between nodes and only sampled once.
    * @param graph Graph to be sampled.
    * @param variant_high If true, use high value for infectiousness of variant.
    * @return Graph with nodes and edges from the input graph sampled.
    */
Graph<Model, MigrationParameters> draw_sample(Graph<Model, MigrationParameters>& graph, bool variant_high);

} // namespace osecirvvs
} // namespace mio

#endif // ODESECIRVVS_PARAMETER_SPACE_H
