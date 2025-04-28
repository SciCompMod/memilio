/* 
* Copyright (C) 2020-2025 MEmilio
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

#ifndef MIO_UTILS_PROFILER_H
#define MIO_UTILS_PROFILER_H

#include "memilio/config.h"

#ifdef MEMILIO_ENABLE_PROFILING

#include "gperftools/profiler.h"

/**
* Start recording a runtime profile from this line until MEMILIO_PROFILER_STOP().
* If using these macros, the `CPUPROFILE` environment variable should not be set,
* the path to the profile output file specified by the macro argument is used instead.
* @param file path to the file where the profile should be stored, as `const char*`.
*/
#define MEMILIO_PROFILER_START(file) ProfilerStart(file)

/**
* Stop recording the runtime profile.
*/
#define MEMILIO_PROFILER_STOP() ProfilerStop()

#else //MEMILIO_ENABLE_PROFILING

//empty definitions so the macros can remain in the code if profiling is disabled.
#define MEMILIO_PROFILER_START(file)
#define MEMILIO_PROFILER_STOP()

#endif //MEMILIO_ENABLE_PROFILING

#endif //MIO_UTILS_PROFILER_H
