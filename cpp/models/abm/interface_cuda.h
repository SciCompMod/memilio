#ifndef MIO_ABM_INTERFACE_CUDA_H
#define MIO_ABM_INTERFACE_CUDA_H

#include <cstdint>
#include <vector>
#include "abm/location_type.h"
#include "abm/infection_state.h"

/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, Khoa Nguyen, Sascha Korf, Carlotta Gerstein
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

namespace mio
{
namespace abm
{

struct GPurson{

        GPurson(LocationType c_loc, uint32_t pid, InfectionState state) : current_loc(c_loc), id(pid), infection_state(state)
        {
            
        }

        LocationType current_loc;
        uint32_t id;
        InfectionState infection_state;

};

// Simple CUDA-compatible structure to hold person data
struct CudaPerson {
    uint32_t id = 0;
    double time_at_location_hours = 0.0;

    CudaPerson(uint32_t p_id, double t)
        : id(p_id), time_at_location_hours(t) {}
};
std::vector<LocationType> mobility_rules(const std::vector<GPurson>& gPursons, int num_persons);
std::vector<double> logTimeAtLocationCuda(const std::vector<CudaPerson>& cuda_persons, int num_persons);
} // namespace abm
} // namespace mio

#endif
