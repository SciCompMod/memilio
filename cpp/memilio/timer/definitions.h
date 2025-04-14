#ifndef MIO_TIMER_DEFINITIONS_H
#define MIO_TIMER_DEFINITIONS_H

#include "memilio/utils/mioomp.h"

#include <chrono>
#include <cstddef>

// clang-format off
#define CONST_LITERAL(literal) []() constexpr { return literal; }
// clang-format on

namespace mio
{
namespace details
{

constexpr size_t constexpr_strlen(const char* string_literal)
{
    // returns length of string, does not count terminating 0 char
    return string_literal[0] == '\0' ? 0 : constexpr_strlen(string_literal + 1) + 1;
}

} // namespace details

namespace timing
{

#ifdef MEMILIO_ENABLE_OPENMP

using TimeType     = decltype(omp_get_wtime());
using DurationType = decltype(omp_get_wtime());

inline double time_in_seconds(DurationType time)
{
    return time;
}

inline TimeType get_time_now()
{
    return omp_get_wtime();
}

#else

using TimeType     = std::chrono::steady_clock::time_point;
using DurationType = std::chrono::steady_clock::duration;

inline double time_in_seconds(DurationType time)
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(time).count();
}

inline TimeType get_time_now()
{
    return std::chrono::steady_clock::now();
}

#endif

} // namespace timing
} // namespace mio

#endif // MIO_TIMER_DEFINITIONS_H
