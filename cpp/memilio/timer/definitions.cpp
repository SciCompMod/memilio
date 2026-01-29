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
#include "memilio/timer/definitions.h"

namespace mio
{
namespace timing
{

#ifdef MEMILIO_ENABLE_OPENMP

double time_in_seconds(DurationType duration)
{
    return duration;
}

TimeType get_time_now()
{
    return omp_get_wtime();
}

#else

double time_in_seconds(DurationType duration)
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(duration).count();
}

TimeType get_time_now()
{
    return std::chrono::steady_clock::now();
}

#endif

} // namespace timing
} // namespace mio
