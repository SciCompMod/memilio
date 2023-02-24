#include "memilio/utils/miompi.h"
#include <mpi.h>

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
