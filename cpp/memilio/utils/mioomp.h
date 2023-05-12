#ifndef MIO_UTIL_OPENMP_H
#define MIO_UTIL_OPENMP_H

#include "memilio/config.h"
#include "memilio/utils/compiler_diagnostics.h"

#ifdef MEMILIO_ENABLE_OPENMP
#include "omp.h"

#define PRAGMA_OMP(x) _Pragma(QUOTE(omp x))

#else

#define PRAGMA_OMP(x)

#endif

#endif
