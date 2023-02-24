
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
#ifdef MEMILIO_ENABLE_MPI
using Comm = MPI_Comm;
#else
using Comm = void*
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
