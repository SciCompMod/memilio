#ifndef MIO_UTILS_PROFILER_H
#define MIO_UTILS_PROFILER_H

#include "memilio/config.h"

#ifdef MEMILIO_ENABLE_PROFILING

#include "gperftools/profiler.h"
#define MEMILIO_PROFILER_START(file) ProfilerStart(file)
#define MEMILIO_PROFILER_STOP() ProfilerStop()

#else //MEMILIO_ENABLE_PROFILING

#define MEMILIO_PROFILER_START(file)
#define MEMILIO_PROFILER_STOP()

#endif //MEMILIO_ENABLE_PROFILING

#endif //MIO_UTILS_PROFILER_H
