/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Martin J. Kuehn
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
