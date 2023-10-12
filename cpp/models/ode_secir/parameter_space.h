/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#ifndef ODESECIR_PARAMETER_SPACE_H
#define ODESECIR_PARAMETER_SPACE_H

#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/memory.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/parameter_distributions.h"
#include "ode_secir/model.h"

#include <assert.h>
#include <string>
#include <vector>
#include <random>
#include <memory>

namespace mio
{
namespace osecir
{
/* Sets alls SecirParams parameters normally distributed, 
*  using the current value and a given standard deviation
* @param[inout] params SecirParams including contact patterns for alle age groups
* @param[in] t0 start time
* @param[in] tmax end time
* @param[in] dev_rel maximum relative deviation from particular value(s) given in params
*/
void set_params_distributions_normal(Model& model, double t0, double tmax, double dev_rel);

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
* as SecirParams parameter values (cf. UncertainValue and SecirParams classes)
* @param[inout] model Model including contact patterns for alle age groups
*/
void draw_sample(Model& model);

Graph<Model, MigrationParameters> draw_sample(Graph<Model, MigrationParameters>& graph);

} // namespace osecir
} // namespace mio

#endif // ODESECIR_PARAMETER_SPACE_H
