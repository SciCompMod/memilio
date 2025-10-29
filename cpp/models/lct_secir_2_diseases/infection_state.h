/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Annika Jungklaus, Lena Ploetzke
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
    // Notation: State_[Infection number][disease]
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
    // Recovered from both diseases
    Recovered_ab = 25,
    Count        = 26
};

} // namespace lsecir2d
} // namespace mio

#endif // LCT_SECIR_2_DISEASES_INFECTIONSTATE_H