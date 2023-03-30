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
#ifndef MIO_UTILS_MPI_H
#define MIO_UTILS_MPI_H

#include "memilio/config.h"

#ifdef MEMILIO_ENABLE_MPI
#include "mpi.h"
#endif

namespace mio
{
namespace mpi
{
/**
* Alias for MPI_Comm, to be used in APIs to avoid
* having to use #ifdef everywhere.
*/
#ifdef MEMILIO_ENABLE_MPI
using Comm = MPI_Comm;
#else
using Comm = void*;
#endif

/**
* Get the global MPI communicator.
*/
Comm get_world();

/**
* Initialize MPI.
*/
void init();

/**
* Finalize MPI.
*/
void finalize();

/**
* Returns true if the calling process is the root process.
*/
bool is_root();

} // namespace mpi
} // namespace mio

#endif
