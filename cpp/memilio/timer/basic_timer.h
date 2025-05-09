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
#ifndef MIO_TIMER_BASIC_TIMER_H
#define MIO_TIMER_BASIC_TIMER_H

#include "memilio/timer/definitions.h"
#include "memilio/utils/logging.h"

namespace mio
{
namespace timing
{

/// @brief A minimal timer class.
class BasicTimer
{
public:
    /// @brief Start the timer. Must be followed by exactly one stop.
    void start()
    {
        should_be_running(false);
        set_running(true);
        m_start_time = get_time_now();
    }

    /// @brief Stop the timer and update the elapsed time. After calling stop, the timer may be started again.
    void stop()
    {
        should_be_running(true);
        set_running(false);
        const TimeType stop_time = get_time_now();
        m_elapsed_time += stop_time - m_start_time;
    }

    /// @brief Set the elapsed time to 0. Only call while the timer is stopped.
    void reset()
    {
        should_be_running(false);
        m_elapsed_time = mio::timing::DurationType{0};
    }

    /// @brief Get the total time spent between starts and stops. Only call while the timer is stopped.
    DurationType get_elapsed_time() const
    {
        should_be_running(false);
        return m_elapsed_time;
    }

    ~BasicTimer()
    {
        should_be_running(false);
    }

private:
    TimeType m_start_time; ///< The last time point at which start() was called
    DurationType m_elapsed_time{0}; ///< The total time spent between starts and stops.

#ifndef NDEBUG
    bool m_is_running = false; ///< In Debug builds, tracks whether the timer is running.

    /// @brief In Debug builds, change whether the timer is running or not.
    void set_running(bool new_state)
    {
        m_is_running = new_state;
    }

    /// @brief In Debug builds, check that the state of m_is_running is as expected. Log an error if not.
    void should_be_running(bool expected) const
    {
        if (m_is_running != expected) {
            mio::log_error("A BasicTimer was {}running while expected to be {}. "
                           "This may be due to mixing up start and stop, calling get_elapsed_time while the timer "
                           "is running, or using the same timer with multiple threads. "
                           "Consider using an AutoTimer with name (and scope) to avoid this.",
                           m_is_running ? "" : "not ", expected ? "started" : "stopped");
        }
        // else: everything ok.
    }
#else
    void set_running(bool)
    {
        // noop
    }
    void should_be_running(bool) const
    {
        // noop
    }
#endif
};

} // namespace timing
} // namespace mio

#endif // MIO_TIMER_BASIC_TIMER_H
