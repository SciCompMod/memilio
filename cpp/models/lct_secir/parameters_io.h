/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Lena Ploetzke
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

#ifndef LCTSECIR_PARAMETERS_IO_H
#define LCTSECIR_PARAMETERS_IO_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "lct_secir/parameters.h"
#include "lct_secir/infection_state.h"
#include "memilio/math/eigen.h"

#include <string>

namespace mio
{
namespace lsecir
{

/**
* @brief Computes an initialization vector for a LCT model with data from RKI.
*   
* Computes an initial value vector for an LCT model with the defined infectionState and the parameters. 
* For the computation expected sojourntime in the subcompartments are used. To calculate the initial values, 
* we assume for simplicity that individuals stay in the subcompartment for exactly the expected time.
* The RKI data are linearly interpolated within one day.
* The RKI data should contain data for each needed day without division of age groups, the completeness of the dates is not verified.
* Data could be downloaded eg with the file pycode/memilio-epidata/memilio/epidata/getCaseData.py, 
* which creates a file named cases_all_germany.json or a similar name. One should set impute_dates=True so that missing dates are imputed.
*
* @param path Path to the RKI file.
* @param date Date for which the initial values should be computed. date is the start date of the simulation.
* @param infectionState InfectionState used to calculate the initial values. Defines eg the size of the calculated vector.
* @param parameters Parameters used to calculate the initial values. Defines eg the expected sojourntimes.
* @param total_population Total size of the population of the considered region. 
*       The sum of all values in the returned vector will be equal to that value.
* @param scale_confirmed_cases Factor by which to scale the confirmed cases of rki data to consider unreported cases.
* @returns Initialization Vector or any io errors that happen during reading of the files.
*/
IOResult<Eigen::VectorXd> get_initial_data_from_file(std::string const& path, Date date, InfectionState infectionState,
                                                     Parameters&& parameters, ScalarType total_population,
                                                     ScalarType scale_confirmed_cases = 1.);

} // namespace lsecir
} // namespace mio
#endif // MEMILIO_HAS_JSONCPP

#endif // LCTSECIR_PARAMETERS_IO_H