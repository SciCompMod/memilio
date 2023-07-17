/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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

#ifndef IDE_SEIR_INFECTIONSTATE_H
#define IDE_SEIR_INFECTIONSTATE_H

namespace mio
{
namespace iseir
{
    /**
    * @brief The InfectionState enum describes the possible categories for the infectious state of persons 
    *    in ide_seir model.
    */
    enum class InfectionState
    {
        S = 0,
        E,
        I,
        R,
        Count = 4
    };

} // namespace iseir
} // namespace mio

#endif
