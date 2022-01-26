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
namespace mio
{
namespace vaccinated
{

    /**
 * @brief The InfectionState enum describes the possible
 * categories for the infectious state of persons
 */
    enum class InfectionState
    {
        Susceptible    = 0,
        SusceptibleV1  = 1,
        Exposed        = 2,
        ExposedV1      = 3,
        ExposedV2      = 4,
        Carrier        = 5,
        CarrierV1      = 6,
        CarrierV2      = 7,
        CarrierT       = 8,
        CarrierTV1     = 9,
        CarrierTV2     = 10,
        Infected       = 11,
        InfectedV1     = 12,
        InfectedV2     = 13,
        InfectedT      = 14,
        InfectedTV1    = 15,
        InfectedTV2    = 16,
        Hospitalized   = 17,
        HospitalizedV1 = 18,
        HospitalizedV2 = 19,
        ICU            = 20,
        ICUV1          = 21,
        ICUV2          = 22,
        Recovered      = 23,
        Dead           = 24,
        InfTotal       = 25,
        Count          = 26
    };

} // namespace vaccinated
} // namespace mio

#endif
