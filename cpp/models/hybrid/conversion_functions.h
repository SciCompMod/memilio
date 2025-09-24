/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Julia Bicker
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

#ifndef MIO_HYBRID_CONVERSION_FUNCTIONS_H
#define MIO_HYBRID_CONVERSION_FUNCTIONS_H

#include "hybrid/temporal_hybrid_model.h"
#include "d_abm/simulation.h"
#include "d_abm/single_well.h"
#include "smm/simulation.h"
#include "memilio/compartments/simulation.h"
#include "ode_secir/model.h"

namespace mio
{
namespace hybrid
{

// This header contains template specilizations for the convert_model function, see mio::hybrid::TemporalHybridSimulation. This function is needed to convert one model to another when the switching condition in the temporal-hybrid model is fulfilled. The TemporalHybridSimulation can be used with any combination of two models, but the template specilizations of the convert_model function have to be provided here.

template <>
void convert_model(const dabm::Simulation<SingleWell<mio::osecir::InfectionState>>& current_model,
                   smm::Simulation<ScalarType, 1, mio::osecir::InfectionState>& target_model);

template <>
void convert_model(const smm::Simulation<ScalarType, 1, mio::osecir::InfectionState>& current_model,
                   dabm::Simulation<SingleWell<mio::osecir::InfectionState>>& target_model);

template <>
void convert_model(const dabm::Simulation<SingleWell<mio::osecir::InfectionState>>& current_model,
                   mio::Simulation<ScalarType, mio::osecir::Model<ScalarType>>& target_model);

template <>
void convert_model(const mio::Simulation<ScalarType, mio::osecir::Model<ScalarType>>& current_model,
                   dabm::Simulation<SingleWell<mio::osecir::InfectionState>>& target_model);

} //namespace hybrid

} //namespace mio

#endif //MIO_HYBRID_CONVERSION_FUNCTIONS_H
