/*
* Copyright (C) 2020-2025 MEmilio
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
#ifndef MIO_TEST_UTILS_H
#define MIO_TEST_UTILS_H

#include "memilio/config.h"
#include "memilio/utils/logging.h"

#include "gtest/gtest.h"
#include "spdlog/sinks/ostream_sink.h"

#include <sstream>

namespace mio
{

/// @brief Can be used to redirect an spdlog::logger. This may cause unintended side effects, like silencing errors!
class RedirectLogger
{
public:
    /**
     * @brief Create a logger that can temporarily capture another's output into a readable/viewable string stream.
     * @param[in] level The LogLevel at which the redirect logger will record logs. Uses "warn" by default.
     */
    RedirectLogger(LogLevel level = LogLevel::warn)
        : m_is_captured(false)
        , m_output()
        , m_sink(std::make_shared<spdlog::sinks::ostream_sink_st>(m_output))
        , m_logger("redirect", m_sink)
        , m_target(nullptr)
    {
        m_logger.set_level(details::get_spdlog_level(level));
    }

    ~RedirectLogger()
    {
        if (m_is_captured) {
            release();
        }
    }

    /**
     * @brief Redirect the output of a logger until it is released.
     * While captured, the logger will write to an internal stream that can be viewed or (destructively) read.
     */
    void capture(std::shared_ptr<spdlog::logger> target = spdlog::default_logger())
    {
        assert(!m_is_captured);
        m_target = target;
        std::swap(*target, m_logger);
        m_is_captured = true;
    }

    /**
     * @brief Restore the logger.
     */
    void release()
    {
        assert(m_is_captured);
        assert(m_target != nullptr);
        std::swap(m_logger, *m_target);
        m_is_captured = false;
    }

    /**
     * @brief Look at the logger content non-destructively.
     */
    const std::string_view view()
    {
        assert(m_is_captured);
        return m_output.view();
    }

    /**
     * @brief Clear the logger content and return it as a string.
     */
    std::string read()
    {
        assert(m_is_captured);
        m_logger.flush();
        std::ostringstream out;
        std::swap(out, m_output);
        return out.str();
    }

private:
    bool m_is_captured; ///< Whether capture is active.
    std::ostringstream m_output; ///< Stream the logger will be redirected to.
    std::shared_ptr<spdlog::sinks::ostream_sink_st> m_sink; ///< Sink linking the stream to the logger.
    spdlog::logger m_logger; ///< The replacement logger. Stores the replaced logger during capture.
    std::shared_ptr<spdlog::logger> m_target; ///< Copy of the target pointer to restore the captured logger.
};

/**
 * @brief Configures the mode of death tests based on the OpenMP configuration
 */
inline void set_death_test_mode()
{
#ifdef MEMILIO_ENABLE_OPENMP
    // If OpenMP is enabled, use "threadsafe" mode to silence gtest warnings.
    GTEST_FLAG_SET(death_test_style, "threadsafe");
#else
    // Without OpenMP, use "fast" mode to avoid GCov / sanitizer related errors.
    GTEST_FLAG_SET(death_test_style, "fast");
#endif
}

} // namespace mio
#endif //MIO_TEST_UTILS_H
