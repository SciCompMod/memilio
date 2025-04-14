#ifndef MEMILIO_TIMER_LIST_PRINTER_H
#define MEMILIO_TIMER_LIST_PRINTER_H

#include "memilio/timer/registration.h"

#include <ios>
#include <iostream>
#include <list>
#include <map>
#include <string>

namespace mio
{
namespace timing
{

struct ListPrinter : public Printer {
    void print(const std::list<TimerRegistration>& timer_register, std::ostream& out) override
    {
        bool is_multithreaded = false;
        const auto indent     = "  ";
        // Write out all timers
        out << "All Timers: " << timer_register.size() << "\n";
        for (const auto& [name, timer, thread] : timer_register) {
            out << indent << name << ": " << std::scientific << time_in_seconds(timer.get_elapsed_time()) << " ("
                << thread << ")\n";
            is_multithreaded |= thread > 0;
        }
        // Write out timers accumulated over threads by name
        if (is_multithreaded) {
            // dedupe list entries from parallel execution
            std::map<std::string, DurationType> deduper;
            for (const auto& [name, timer, thread] : timer_register) {
                deduper[name] += timer.get_elapsed_time();
            }
            out << "Unique Timers (accumulated): " << deduper.size() << "\n";
            for (const auto& [name, time] : deduper) {
                out << indent << name << ": " << std::scientific << time_in_seconds(time) << "\n";
            }
        }
    }
};

} // namespace timing
} // namespace mio

#endif // MEMILIO_TIMER_LIST_PRINTER_H
