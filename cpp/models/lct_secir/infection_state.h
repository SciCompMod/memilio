/* 
* Copyright (C) 2020-2025 MEmilio
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

#ifndef LCT_SECIR_INFECTIONSTATE_H
#define LCT_SECIR_INFECTIONSTATE_H

namespace mio
{
namespace lsecir
{

/**
 * @brief The InfectionState enum describes the basic
 * categories for the infection state of persons.
 */
enum class InfectionState
{
    Susceptible        = 0,
    Exposed            = 1,
    InfectedNoSymptoms = 2,
    InfectedSymptoms   = 3,
    InfectedSevere     = 4,
    InfectedCritical   = 5,
    Recovered          = 6,
    Dead               = 7,
    Count              = 8
};

/**
 * @brief The InfectionTransition enum describes the possible
 * transitions of the infectious state of persons.
 */
enum class InfectionTransition
{
    SusceptibleToExposed                 = 0,
    ExposedToInfectedNoSymptoms          = 1,
    InfectedNoSymptomsToInfectedSymptoms = 2,
    InfectedNoSymptomsToRecovered        = 3,
    InfectedSymptomsToInfectedSevere     = 4,
    InfectedSymptomsToRecovered          = 5,
    InfectedSevereToInfectedCritical     = 6,
    InfectedSevereToRecovered            = 7,
    InfectedCriticalToDead               = 8,
    InfectedCriticalToRecovered          = 9,
    Count                                = 10
};

} // namespace lsecir
} // namespace mio

#endif // LCT_SECIR_INFECTIONSTATE_H