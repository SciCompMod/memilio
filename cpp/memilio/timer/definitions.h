/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding
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
#ifndef MIO_TIMER_DEFINITIONS_H
#define MIO_TIMER_DEFINITIONS_H

#include "memilio/utils/mioomp.h"

#include <chrono>
#include <cstddef>

namespace mio
{

namespace timing
{

#ifdef MEMILIO_ENABLE_OPENMP

using TimeType     = decltype(omp_get_wtime());
using DurationType = decltype(omp_get_wtime());

/**
 * @brief Convert a duration to a (floating point) number of seconds.
 * @param[in] duration Any DurationType value, mainly `BasicTimer::get_elapsed_time()`.  
 * @return The duration in seconds.
 */
inline double time_in_seconds(DurationType duration)
{
    return duration;
}

/**
 * @brief Get the current time.
 * @return Returns omp_get_wtime() if OpenMP is enabled, stead_clock::now() otherwise.
 */
inline TimeType get_time_now()
{
    return omp_get_wtime();
}

#else

using TimeType     = std::chrono::steady_clock::time_point;
using DurationType = std::chrono::steady_clock::duration;

/**
 * @brief Convert a duration to a (floating point) number of seconds.
 * @param[in] duration Any DurationType value, mainly `BasicTimer::get_elapsed_time()`.  
 * @return The duration in seconds.
 */
inline double time_in_seconds(DurationType duration)
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(duration).count();
}

/**
 * @brief Get the current time.
 * @return Returns omp_get_wtime() if OpenMP is enabled, stead_clock::now() otherwise.
 */
inline TimeType get_time_now()
{
    return std::chrono::steady_clock::now();
}

#endif

} // namespace timing
} // namespace mio

#endif // MIO_TIMER_DEFINITIONS_H
