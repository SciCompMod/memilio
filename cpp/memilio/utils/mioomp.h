/* 
* Copyright (C) 2020-2026 MEmilio
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

#include "memilio/config.h" // included for MEMILIO_ENABLE_OPENMP
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
namespace omp
{

/// @brief Get an id for the current OpenMP thread. When OpenMP is disabled, this is always 0.
int get_thread_id();

/// @brief Get the current number of OpenMP threads. When OpenMP is disabled, this is always 1.
int get_num_threads();

/// @brief Get an upper bound to the number of OpenMP threads. When OpenMP is disabled, this is always 1.
int get_max_threads();

} // namespace omp
} // namespace mio

#endif
