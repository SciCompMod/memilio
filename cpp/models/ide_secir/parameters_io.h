/*
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Lena Ploetzke, Anna Wendler
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

#include "ide_secir/model_ide.h"
#include "memilio/io/io.h"
#include "memilio/utils/date.h"

#include <string>

namespace mio
{
namespace isecir
{

/**
* @brief Computes a TimeSeries of flows to provide initial data for an IDE-SECIR model with data from RKI.
*   
* The flows InfectedNoSymptomsToInfectedSymptoms are calculated using the confirmed cases in the RKI data.
* If necessary, the RKI data are linearly interpolated within one day.
* The RKI data should contain data for each needed day without division of age groups, the completeness of the dates is not verified.
* Data can be downloaded e.g. with the file pycode/memilio-epidata/memilio/epidata/getCaseData.py, 
* which creates a file named cases_all_germany.json or a similar name. One should set impute_dates=True so that missing dates are imputed.
*
* The flows InfectedSymptomsToInfectedSevere, InfectedSymptomsToRecovered, InfectedSevereToInfectedCritical,
* InfectedSevereToRecovered, InfectedCriticalToDead and InfectedCriticalToRecovered can then be calculated using
* the InfectedNoSymptomsToInfectedSymptoms flow with the standard formula from the IDE model.
* The ExposedToInfectedNoSymptoms and InfectedNoSymptomsToInfectedSymptoms flows are calculated 
* using the means of the respective TransitionDistribution. 
* The flow InfectedNoSymptomsToInfectedSymptoms is calculated with the standard formula from the IDE model
* using the results for ExposedToInfectedNoSymptoms.
*
* The number of deaths used in the model is set to the number given in the RKI data.
* We also set the number of total confirmed cases in the model. 
* Therefore the initialization method using the total confirmed cases is used in the model. See also the documentation of the model class.
* 
* The start date of the model simulation is set to t0=0.
*
* @param[in, out] model The model for which the initial flows should be computed.
* @param[in] dt Time step size.
* @param[in] path Path to the RKI file.
* @param[in] date The start date of the simulation and the last time point of the TimeSeries used for initialization.
* @param[in] scale_confirmed_cases Factor by which to scale the confirmed cases of rki data to consider unreported cases.
* @returns Any io errors that happen during reading of the files.
*/
IOResult<void> set_initial_flows(Model& model, ScalarType dt, std::string const& path, Date date,
                                 ScalarType scale_confirmed_cases = 1.);

} // namespace isecir
} // namespace mio

#endif // MEMILIO_HAS_JSONCPP

#endif // IDE_INITIALFLOWS_H