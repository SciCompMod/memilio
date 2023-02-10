/* 
* Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
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
#ifndef MIO_EPI_SIMULATION_DAY_H
#define MIO_EPI_SIMULATION_DAY_H

#include "memilio/utils/index.h"

namespace mio
{

/**
* Represents the simulation time as an integer index.
*/
class SimulationDay : public Index<SimulationDay>
{
public:
    using Index<SimulationDay>::Index;
};

} // namespace mio

#endif //MIO_EPI_SIMULATION_DAY_H
