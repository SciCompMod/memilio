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
#ifndef MIO_TIMER_LIST_PRINTER_H
#define MIO_TIMER_LIST_PRINTER_H

#include "memilio/timer/registration.h"
#include "memilio/utils/compiler_diagnostics.h"

#include <ostream>
#include <list>
#include <map>
#include <string>

namespace mio
{
namespace timing
{

struct ListPrinter : public Printer {
    /**
     * @brief Print the elapsed time of all registered timers in a list.
     * If multiple threads were used for timing, a second list with the results aggregated over threads is printed.
     * @param[in] timer_register TimerRegistrar's register of timers.
     * @param[inout] out The stream to write to.
     */
    void print(const std::list<TimerRegistration>& timer_register, std::ostream& out) override
    {
        bool is_multithreaded = false;
        const auto indent     = "  ";
        // Write out all timers
        out << "All Timers: " << timer_register.size() << "\n";
        for (const auto& [name, scope, timer, thread, rank] : timer_register) {
            out << indent << qualified_name(name, scope) << ": " << std::scientific
                << time_in_seconds(timer.get_elapsed_time()) << " (" << rank << ", " << thread << ")\n";
            is_multithreaded |= thread > 0 || rank > 0;
        }
        // Write out timers accumulated over threads by name
        if (is_multithreaded) {
            // dedupe list entries from parallel execution
            std::map<std::string, DurationType> deduper;
            for (const auto& [name, scope, timer, thread, rank] : timer_register) {
                deduper[qualified_name(name, scope)] += timer.get_elapsed_time();
            }
            out << "Unique Timers (accumulated): " << deduper.size() << "\n";
            for (const auto& [q_name, time] : deduper) {
                out << indent << q_name << ": " << std::scientific << time_in_seconds(time) << "\n";
            }
        }
    }
};

} // namespace timing
} // namespace mio

#endif // MIO_TIMER_LIST_PRINTER_H
