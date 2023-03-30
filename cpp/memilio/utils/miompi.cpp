/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#include "memilio/utils/miompi.h"

#ifdef MEMILIO_ENABLE_MPI
#include <mpi.h>
#endif

namespace mio
{
namespace mpi
{
Comm get_world()
{
#ifdef MEMILIO_ENABLE_MPI
    static Comm world = MPI_COMM_WORLD;
#else
    static Comm world = nullptr;
#endif
    return world;
}

void init()
{
#ifdef MEMILIO_ENABLE_MPI
    MPI_Init(nullptr, nullptr);
#endif
}

void finalize()
{
#ifdef MEMILIO_ENABLE_MPI
    MPI_Finalize();
#endif
}

bool is_root()
{
#ifdef MEMILIO_ENABLE_MPI
    int rank;
    MPI_Comm_rank(get_world(), &rank);
    return rank == 0;
#endif
    return true;
}

} // namespace mpi
} // namespace mio
