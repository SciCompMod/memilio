#ifndef MIO_UTIL_OPENMP_H
#define MIO_UTIL_OPENMP_H

#include "memilio/config.h"
#include "memilio/utils/compiler_diagnostics.h"

#ifdef MEMILIO_ENABLE_OPENMP
#include "omp.h"

/**
* Macro for OpenMP directives.
* OpenMP is enabled.
* `PRAGMA_OMP(parallel for)` evaluates to `#pragma omp parallel for`
*/
#define PRAGMA_OMP(x) _Pragma(QUOTE(omp x))

#else

/**
* Macro for OpenMP directives. 
* Evaluates to nothing because OpenMP is disabled.
* Unknown pragmas are ignored by the compiler, but the compiler emits warnings, so
* we remove the #pragma completely.
*/
#define PRAGMA_OMP(x)

#endif

namespace mio
{

inline int get_omp_thread_id()
{
#ifdef MEMILIO_ENABLE_OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

inline int get_omp_num_threads()
{
#ifdef MEMILIO_ENABLE_OPENMP
    return omp_get_num_threads();
#else
    return 1;
#endif
}

inline int get_omp_max_threads()
{
#ifdef MEMILIO_ENABLE_OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

} // namespace mio

#endif
