/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#ifndef ODE_MSEIRS4_INFECTION_STATE_H
#define ODE_MSEIRS4_INFECTION_STATE_H

namespace mio
{
namespace omseirs4
{
enum class InfectionState
{
    MaternalImmune = 0,
    S1,
    S2,
    S3,
    S4,
    E1,
    E2,
    E3,
    E4,
    I1,
    I2,
    I3,
    I4,
    R1,
    R2,
    R3,
    R4,
    Count
};

} // namespace omseirs4
} // namespace mio

#endif // ODE_MSEIRS4_INFECTION_STATE_H
