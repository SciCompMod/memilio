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
#ifndef MIO_TIMER_AUTO_TIMER_H
#define MIO_TIMER_AUTO_TIMER_H

#include "memilio/timer/basic_timer.h"
#include "memilio/timer/named_timer.h"
#include "memilio/utils/string_literal.h"

namespace mio
{
namespace timing
{

/**
 * @brief Timer that automatically starts when it is created, and stops when it is destroyed.
 * @tparam Name, Scope The name and scope of a NamedTimer. Do not set these if you want to use a BasicTimer.
 */
template <StringLiteral Name, StringLiteral Scope = "">
class AutoTimer
{
public:
    /// @brief Run the NamedTimer given by the template parameter(s) Name (and Scope).
    AutoTimer()
        : m_timer(NamedTimer<Name, Scope>::get_instance())
    {
        m_timer.start();
    }

    /// @brief Run the given BasicTimer. Does not take ownership, so mind the timer's lifetime!
    AutoTimer(BasicTimer& timer)
        : m_timer(timer)
    {
        static_assert(Name.empty() && Scope.empty(),
                      "Do not set the Name and Scope templates when using this constructor.");
        m_timer.start();
    }

    AutoTimer(AutoTimer&)  = delete;
    AutoTimer(AutoTimer&&) = delete;

    ~AutoTimer()
    {
        m_timer.stop();
    }

private:
    BasicTimer& m_timer; ///< Reference to the timer so it can be stopped in AutoTimer's destructor.
};

// Deduction guide that allows omitting the template parameter when using the BasicTimer constructor.
AutoTimer(BasicTimer& timer) -> AutoTimer<"">;

} // namespace timing
} // namespace mio

#endif // MIO_TIMER_AUTO_TIMER_H
