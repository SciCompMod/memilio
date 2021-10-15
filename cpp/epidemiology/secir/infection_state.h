/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#ifndef INFECTIONSTATE_H
#define INFECTIONSTATE_H

namespace epi
{

/**
 * @brief The InfectionState enum describes the possible
 * categories for the infectious state of persons
 */
enum class InfectionState
{
    Susceptible  = 0,
    Exposed      = 1,
    Carrier      = 2,
    Infected     = 3,
    Hospitalized = 4,
    ICU          = 5,
    Recovered    = 6,
    Dead         = 7,
    Count        = 8
};

/**
 * @brief The InfectionState enum describes the possible
 * categories for the infectious state of persons
 */
enum class InfectionStateV
{
    Susceptible   = 0,
    SusceptibleV1 = 1,
    //SusceptibleV2  = 2,
    Exposed   = 2,
    ExposedV1 = 3,
    //ExposedV2      = 5,
    Carrier   = 4,
    CarrierV1 = 5,
    //CarrierV2      = 8,
    CarrierT   = 6,
    CarrierTV1 = 7,
    //CarrierTV2     = 11,
    Infected   = 8,
    InfectedV1 = 9,
    //InfectedV2     = 14,
    InfectedT   = 10,
    InfectedTV1 = 11,
    //InfectedTV2    = 17,
    Hospitalized   = 12,
    HospitalizedV1 = 13,
    //HospitalizedV2 = 20,
    ICU   = 14,
    ICUV1 = 15,
    //ICUV2          = 23,
    Recovered = 16,
    Dead      = 17,
    Count     = 18
};
} // namespace epi

#endif
