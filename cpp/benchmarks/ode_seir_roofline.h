/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Henrik Zunker
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
#ifndef MIO_ODE_SEIR_ROOFLINE_H
#define MIO_ODE_SEIR_ROOFLINE_H

#include <atomic>
#include <cstdlib>
#include <initializer_list>
#include <stdexcept>
#include <string>
#ifdef MEMILIO_ROOFLINE_LIKWID
#include <likwid.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#endif

namespace mio::benchmark_roofline
{
// This opt-in is independent of the historical LIKWID kernel markers. Normal
// builds contain neither counter calls nor GPU profiling APIs in their hot path.
inline bool enabled()
{
    const char* value = std::getenv("MEMILIO_ROOFLINE");
    if (!value || std::string(value) == "0")
        return false;
    if (std::string(value) != "1")
        throw std::invalid_argument("MEMILIO_ROOFLINE must be 0 or 1.");
#ifndef MEMILIO_BENCHMARK_ROOFLINE
    throw std::runtime_error("Roofline requires a separate MEMILIO_BENCHMARK_ROOFLINE=ON build.");
#else
    return true;
#endif
}

inline int positive_setting(const char* key, int fallback)
{
    const char* value = std::getenv(key);
    if (!value)
        return fallback;
    const std::string text(value);
    if (text.empty() || text.size() > 3 || text.find_first_not_of("0123456789") != std::string::npos)
        throw std::invalid_argument(std::string(key) + " must be an integer in 1..366.");
    const int result = std::stoi(text);
    if (result < 1 || result > 366)
        throw std::invalid_argument(std::string(key) + " must be an integer in 1..366.");
    return result;
}

inline int selected_day()
{
    const int day = positive_setting("ROOFLINE_DAY", 2);
    if (day > positive_setting("RUNTIME_DAYS", 32))
        throw std::invalid_argument("ROOFLINE_DAY must not exceed RUNTIME_DAYS.");
    return day; // Human-readable one-based day; the preceding days are unprofiled.
}

#ifdef MEMILIO_BENCHMARK_ROOFLINE
inline std::atomic<bool> cpu_marker_failed{false};

class CpuSession
{
public:
    explicit CpuSession(std::initializer_list<const char*> regions)
        : active(enabled())
    {
        if (!active)
            return;
#ifdef MEMILIO_ROOFLINE_LIKWID
        if (!std::getenv("LIKWID_FILEPATH"))
            throw std::runtime_error("CPU roofline must be launched through likwid-perfctr -m.");
        cpu_marker_failed.store(false);
        likwid_markerInit();
#pragma omp parallel
        {
            likwid_markerThreadInit();
            for (const auto* region : regions)
                if (likwid_markerRegisterRegion(region) != 0)
                    cpu_marker_failed.store(true);
        }
        if (cpu_marker_failed.load()) {
            likwid_markerClose();
            active = false;
            throw std::runtime_error("LIKWID region registration failed; check CPU support and counter access.");
        }
#else
        (void)regions;
        throw std::runtime_error("CPU roofline requires MEMILIO_ROOFLINE_LIKWID=ON.");
#endif
    }
    ~CpuSession()
    {
#ifdef MEMILIO_ROOFLINE_LIKWID
        if (active)
            likwid_markerClose();
#endif
    }
    CpuSession(const CpuSession&) = delete;
    CpuSession& operator=(const CpuSession&) = delete;
    void check() const
    {
        if (active && cpu_marker_failed.load())
            throw std::runtime_error("LIKWID marker collection failed; no valid roofline sample was produced.");
    }
private:
    bool active;
};

template <class Work>
inline void cpu_region(const char* name, Work&& work)
{
#ifdef MEMILIO_ROOFLINE_LIKWID
    // All threads call markers. Align boundaries so a shared memory-controller
    // counter never sees another thread begin the next phase before its stop.
#pragma omp barrier
    if (likwid_markerStartRegion(name) != 0)
        cpu_marker_failed.store(true);
#pragma omp barrier
    work();
#pragma omp barrier
    if (likwid_markerStopRegion(name) != 0)
        cpu_marker_failed.store(true);
#pragma omp barrier
#else
    (void)name;
    (void)work;
    throw std::runtime_error("CPU roofline requires LIKWID support.");
#endif
}
#endif
} // namespace mio::benchmark_roofline
#endif
