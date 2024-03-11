/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Anna Wendler, Lena Ploetzke
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
#ifndef IDE_INITIALFLOWS_H
#define IDE_INITIALFLOWS_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "ide_secir/model.h"
#include "memilio/math/eigen.h"
#include "memilio/io/io.h"
#include "memilio/utils/date.h"

#include <string>

namespace mio
{
namespace isecir
{

ScalarType compute_mean(Eigen::Index idx_CurrentFlow);

void compute_previous_flows(Model& model, Eigen::Index idx_CurrentFlow, Eigen::Index idx_OutgoingFlow,
                            Eigen::Index time_series_index, ScalarType dt);

void compute_flows_with_mean(Model& model, Eigen::Index idx_InfectionTransitions, Eigen::Index idx_OutgoingFlow,
                             ScalarType dt, Eigen::Index current_time_index);

/**
* @brief Computes a TimeSeries of flows to provide initial data for an IDE SECIR model with data from RKI.
*   
* TODO: describe how this is computed
*
* @param[in, out] Model The model for which the inital flows should be computed.
* @param[in] dt Time step size.
* @param[in] path Path to the RKI file.
* @param[in] date The start date of the simulation and the last time point of the TimeSeries used for initialization.
* @param[in] scale_confirmed_cases Factor by which to scale the confirmed cases of rki data to consider unreported cases.
* @returns Any io errors that happen during reading of the files.
*/
IOResult<void> set_initial_flows(Model& model, ScalarType dt, std::string const& path, Date date,
                                 ScalarType scale_confirmed_cases = 1.0);

} // namespace isecir
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // IDE_INITIALFLOWS_H