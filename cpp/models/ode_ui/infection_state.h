/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn
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

#ifndef ODEUI_INFECTIONSTATE_H
#define ODEUI_INFECTIONSTATE_H

namespace mio
{
namespace oui
{

/**
     * @brief The InfectionState enum describes the possible
     * categories for the infectious state of persons
     */
enum class InfectionState
{   
    Susceptible,
    InfectedV1,
    InfectedV2,
    InfectedV3,
    Count
};

} // namespace oui
} // namespace mio

#endif // ODEUI_INFECTIONSTATE_H
