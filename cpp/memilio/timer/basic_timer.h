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
        m_start_time = get_time_now();
    }

    /// @brief Stop the timer and update the elapsed time. After calling stop, the timer may be started again.
    void stop()
    {
        const TimeType stop_time = get_time_now();
        m_elapsed_time += stop_time - m_start_time;
    }

    /// @brief Set the elapsed time to 0. Only call while stopped.
    void reset()
    {
        m_elapsed_time = 0;
    }

    /// @brief Get the total time spent between starts and stops. Do not use while the timer is started.
    DurationType get_elapsed_time() const
    {
        return m_elapsed_time;
    }

private:
    TimeType m_start_time; ///< The last time point at which start() was called
    DurationType m_elapsed_time{0}; ///< The total time spent between starts and stops.
};

} // namespace timing
} // namespace mio

#endif // MIO_TIMER_BASIC_TIMER_H
