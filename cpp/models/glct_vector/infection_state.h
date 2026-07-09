/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Annika Jungklaus
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

#ifndef MIO_GLCT_VECTOR_INFECTIONSTATE_H
#define MIO_GLCT_VECTOR_INFECTIONSTATE_H

namespace mio
{
namespace glvector
{

/// @brief The InfectionState enum describes the basic categories for the infection state of persons and vectors.
enum class InfectionState
{
    UnbornHuman             = 0,
    SusceptibleHuman        = 1,
    ExposedHuman            = 2,
    InfectedHuman            = 3,
    RecoveredHuman          = 4,
    DeadHuman               = 5,
    NaturalDeathHuman      = 6,
    UnbornVector            = 7,
    AquaticVector             = 8,
    SusceptibleVector       = 9,
    ExposedVector          = 10,
    TransmitterVector      = 11,
    DeadVector             = 12,
    Count                   = 13

};

} // namespace glvector
} // namespace mio

#endif // MIO_GLCT_VECTOR_INFECTIONSTATE_H