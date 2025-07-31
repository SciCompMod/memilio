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

#ifndef LCT_SECIR_2_DISEASES_INFECTIONSTATE_H
#define LCT_SECIR_2_DISEASES_INFECTIONSTATE_H

namespace mio
{
namespace lsecir2d
{

/**
 * @brief The InfectionState enum describes the basic
 * categories for the infection state of persons.
 */
enum class InfectionState
{
    Susceptible = 0,
    // State_[Infection number][disease]
    // first infection with disease a
    Exposed_1a            = 1,
    InfectedNoSymptoms_1a = 2,
    InfectedSymptoms_1a   = 3,
    InfectedSevere_1a     = 4,
    InfectedCritical_1a   = 5,
    Recovered_a           = 6,
    Dead_a                = 7,
    // second infection with disease a
    Exposed_2a            = 8,
    InfectedNoSymptoms_2a = 9,
    InfectedSymptoms_2a   = 10,
    InfectedSevere_2a     = 11,
    InfectedCritical_2a   = 12,
    // R and D for disease a

    // first infection with disease b
    Exposed_1b            = 13,
    InfectedNoSymptoms_1b = 14,
    InfectedSymptoms_1b   = 15,
    InfectedSevere_1b     = 16,
    InfectedCritical_1b   = 17,
    Recovered_b           = 18,
    Dead_b                = 19,
    // second infection with disease b
    Exposed_2b            = 20,
    InfectedNoSymptoms_2b = 21,
    InfectedSymptoms_2b   = 22,
    InfectedSevere_2b     = 23,
    InfectedCritical_2b   = 24,
    // R and D for disease b

    // Recovered from both diseases
    Recovered_ab = 25,
    Count        = 26
};

/**
 * @brief The InfectionTransition enum describes the possible
 * transitions of the infectious state of persons.
 */
enum class InfectionTransition
{
    // first infection with a
    SusceptibleToExposed_1a                    = 0,
    Exposed_1aToInfectedNoSymptoms_1a          = 1,
    InfectedNoSymptoms_1aToInfectedSymptoms_1a = 2,
    InfectedNoSymptoms_1aToRecovered_a         = 3,
    InfectedSymptoms_1aToInfectedSevere_1a     = 4,
    InfectedSymptoms_1aToRecovered_a           = 5,
    InfectedSevere_1aToInfectedCritical_1a     = 6,
    InfectedSevere_1aToRecovered_a             = 7,
    InfectedCritical_1aToDead_a                = 8,
    InfectedCritical_1aToRecovered_a           = 9,
    // second infection with a
    Recovered_bToExposed_2a                    = 10,
    Exposed_2aToInfectedNoSymptoms_2a          = 11,
    InfectedNoSymptoms_2aToInfectedSymptoms_2a = 12,
    InfectedNoSymptoms_2aToRecovered_ab        = 13,
    InfectedSymptoms_2aToInfectedSevere_2a     = 14,
    InfectedSymptoms_2aToRecovered_ab          = 15,
    InfectedSevere_2aToInfectedCritical_2a     = 16,
    InfectedSevere_2aToRecovered_ab            = 17,
    InfectedCritical_2aToDead_a                = 18,
    InfectedCritical_2aToRecovered_ab          = 19,
    // first infection with b
    SusceptibleToExposed_1b                    = 20,
    Exposed_1bToInfectedNoSymptoms_1b          = 21,
    InfectedNoSymptoms_1bToInfectedSymptoms_1b = 22,
    InfectedNoSymptoms_1bToRecovered_b         = 23,
    InfectedSymptoms_1bToInfectedSevere_1b     = 24,
    InfectedSymptoms_1bToRecovered_b           = 25,
    InfectedSevere_1bToInfectedCritical_1b     = 26,
    InfectedSevere_1bToRecovered_b             = 27,
    InfectedCritical_1bToDead_b                = 28,
    InfectedCritical_1bToRecovered_b           = 29,
    // second infection with b
    Recovered_aToExposed_2b                    = 30,
    Exposed_2bToInfectedNoSymptoms_2b          = 31,
    InfectedNoSymptoms_2bToInfectedSymptoms_2b = 32,
    InfectedNoSymptoms_2bToRecovered_ab        = 33,
    InfectedSymptoms_2bToInfectedSevere_2b     = 34,
    InfectedSymptoms_2bToRecovered_ab          = 35,
    InfectedSevere_2bToInfectedCritical_2b     = 36,
    InfectedSevere_2bToRecovered_ab            = 37,
    InfectedCritical_2bToDead_b                = 38,
    InfectedCritical_2bToRecovered_ab          = 39,
    Count                                      = 40
};

} // namespace lsecir2d
} // namespace mio

#endif // LCT_SECIR_2_DISEASES_INFECTIONSTATE_H