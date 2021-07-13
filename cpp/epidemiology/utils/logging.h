#ifndef LOGGING_H
#define LOGGING_H

#ifdef NDEBUG
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_INFO
#else
#define SPDLOG_ACTIVE_LEVEL SPDLOG_LEVEL_TRACE
#endif

#include "epidemiology/utils/compiler_diagnostics.h"
#include <spdlog/spdlog.h>

namespace epi
{

enum class LogLevel
{
    trace,
    debug,
    info,
    warn,
    err,
    critical,
    off
};

namespace details
{
    inline spdlog::level::level_enum get_spdlog_level(LogLevel level)
    {
        spdlog::level::level_enum l;
        switch (level) {
        case LogLevel::trace:
            l = spdlog::level::trace;
            break;
        case LogLevel::debug:
            l = spdlog::level::debug;
            break;
        case LogLevel::info:
            l = spdlog::level::info;
            break;
        case LogLevel::warn:
            l = spdlog::level::warn;
            break;
        case LogLevel::err:
            l = spdlog::level::err;
            break;
        case LogLevel::critical:
            l = spdlog::level::critical;
            break;
        case LogLevel::off:
            l = spdlog::level::off;
            break;
        default:
            l = spdlog::level::info;
            assert(false && "Unknown LogLevel.");
        }
        return l;
    }
} // namespace details

/**
 * @brief Sets the verbosity of the logger
 */
inline void set_log_level(LogLevel level)
{
    spdlog::set_level(details::get_spdlog_level(level));
}

template <typename... Args> inline void log_info(spdlog::string_view_t fmt, const Args&... args)
{
    spdlog::default_logger_raw()->info(fmt, args...);
}

template <typename... Args> inline void log_error(spdlog::string_view_t fmt, const Args&... args)
{
    spdlog::default_logger_raw()->error(fmt, args...);
}

template <typename... Args> inline void log_warning(spdlog::string_view_t fmt, const Args&... args)
{
    spdlog::default_logger_raw()->warn(fmt, args...);
}

template <typename... Args> inline void log_debug(spdlog::string_view_t fmt, const Args&... args)
{
#ifndef NDEBUG
    spdlog::default_logger_raw()->debug(fmt, args...);
#else 
    unused(fmt, args...);
#endif
}

template <typename... Args>
inline void log(LogLevel level, spdlog::string_view_t fmt, const Args&... args)
{
    spdlog::default_logger_raw()->log(details::get_spdlog_level(level), fmt, args...);
}

} // namespace epi

#endif // LOGGING_H
