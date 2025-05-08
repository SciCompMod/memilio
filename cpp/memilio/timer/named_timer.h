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
#ifndef MIO_TIMER_NAMED_TIMER_H
#define MIO_TIMER_NAMED_TIMER_H

#include "memilio/utils/mioomp.h"
#include "memilio/timer/basic_timer.h"
#include "memilio/timer/timer_registrar.h"
#include "memilio/utils/string_literal.h"

#include <string>

namespace mio
{
namespace timing
{

/// @brief Thread local singleton timer, identified by its name. Best used via AutoTimer.
template <StringLiteral Name, StringLiteral Scope = "">
class NamedTimer : public BasicTimer
{
public:
    /// @brief Do not allow construction other than through get_instance and the default constructor.
    NamedTimer(NamedTimer&)  = delete;
    NamedTimer(NamedTimer&&) = delete;

    /// @brief Get the timer's name as a runtime string.
    static std::string name()
    {
        return std::string(Name);
    }

    /// @brief Get the timer's scope as a runtime string.
    static std::string scope()
    {
        return std::string(Scope);
    }

    /**
     * @brief Get a thread local instance of this timer. It is automatically registered to the TimerRegistrar.
     * 
     * This is the only way to obtain a NamedTimer, as its constructors are disabled. Since there is exactly one timer
     * for each Name and Scope combination, this timer is inherently threadsafe. The Name and Scope together effectively
     * work as a key for a global map, though the map access is made during compilation.
     */
    static NamedTimer& get_instance()
    {
        // create and return a timer
        // the constructor is only called the first time this function is used (for each thread)
        thread_local static NamedTimer t;
        return t;
    }

private:
    /**
     * @brief Create a NamedTimer and register it.
     * Since this is the only available constructor, and the registration is made before construction is done, it is
     * guaranteed (by standard), that NamedTimer is destroyed after TimerRegistrar. See
     * https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2017/n4659.pdf 15.6.2 and 15.4-9.
     */
    NamedTimer()
    {
        TimerRegistrar::get_instance().add_timer({name(), scope(), *this, mio::get_omp_thread_id()});
    }
};

} // namespace timing
} // namespace mio

#endif // MIO_TIMER_NAMED_TIMER_H
