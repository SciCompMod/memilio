/* 
* Copyright (C) 2020-2026 MEmilio
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
#include <ctime>

namespace mio
{
namespace timing
{

#ifdef MEMILIO_ENABLE_OPENMP

using TimeType     = decltype(omp_get_wtime());
using DurationType = decltype(omp_get_wtime());

#else

using TimeType     = std::chrono::steady_clock::time_point;
using DurationType = std::chrono::steady_clock::duration;

#endif

namespace details
{

/// @brief Convert a duration to integer ticks. Useful for serialization.
inline decltype(auto) convert_to_ticks(DurationType duration)
{
#ifdef MEMILIO_ENABLE_OPENMP
    return duration;
#else
    return duration.count();
#endif
}

} // namespace details

/**
 * @brief Convert a duration to a (floating point) number of seconds.
 * @param[in] duration Any DurationType value, mainly `BasicTimer::get_elapsed_time()`.  
 * @return The duration in seconds.
 */
double time_in_seconds(DurationType duration);

/**
 * @brief Get the current time.
 * @return Returns omp_get_wtime() if OpenMP is enabled, steady_clock::now() otherwise.
 */
TimeType get_time_now();

} // namespace timing
} // namespace mio

#endif // MIO_TIMER_DEFINITIONS_H
