#ifndef MIO_TIMER_TIMER_REGISTRAR_H
#define MIO_TIMER_TIMER_REGISTRAR_H

#include "memilio/timer/definitions.h"
#include "memilio/timer/basic_timer.h"

#include <iostream>
#include <list>
#include <map>
#include <mutex>
#include <string>

namespace mio
{
namespace timing
{

struct TimerRegistration {
    std::string name;
    BasicTimer& timer;
    int thread_id;
};

// holds raw pointers to names and durations of the thread_local Timer instances
class TimerRegistrar
{
public:
    // can either use get_instance scheme, or nake all members static. the latter does not use ctor/dtor
    //
    // get_instance scheme has the benefit of guaranteeing correct destruction order, see
    // https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2017/n4659.pdf 15.6.2 and 15.4-9
    // this is due to get_instance being called withing the Timer ctor body, hence TimerRegistrar _finishes_ its ctor
    // before Timer, even if Timer() is called first.
    static TimerRegistrar& get_instance()
    {
        static TimerRegistrar t;
        return t;
    }

    void add_timer(TimerRegistration&& registration)
    {
        m_registration_lock.lock();
        m_register.emplace_back(registration);
        m_registration_lock.unlock();
    }

    const auto& get_register()
    {
        return m_register;
    }

    // questionable whether this should be allowed to be used over print_on_exit.
    void print_timers(std::ostream& out = std::cout)
    {
        PRAGMA_OMP(single)
        {
            const auto indent = "  ";
            out << "All Timers: " << m_register.size() << "\n";
            for (const auto& entry : m_register) {
                out << indent << entry.name << ": " << std::scientific
                    << time_in_seconds(entry.timer.get_elapsed_time()) << "\n";
            }
            // dedupe list entries from parallel execution
            std::map<std::string, DurationType> deduper;
            for (const auto& entry : m_register) {
                deduper[entry.name] += entry.timer.get_elapsed_time();
            }
            out << "Unique Timers (accumulated): " << deduper.size() << "\n";
            for (const auto& entry : deduper) {
                out << indent << entry.first << ": " << std::scientific << time_in_seconds(entry.second) << "\n";
            }
        }
    }

    void print_on_exit() // requires dtor (i.e. get_instance scheme)
    {
        m_print_on_death = true;
    }

private:
    TimerRegistrar() = default;

    ~TimerRegistrar()
    {
        if (m_print_on_death) {
            PRAGMA_OMP(single)
            {
                std::cout << "Final Timer Summary\n";
                print_timers();
            }
        }
    }

    std::list<TimerRegistration> m_register;
    std::mutex m_registration_lock;
    bool m_print_on_death;
};

} // namespace timing
} // namespace mio

#endif // MIO_TIMER_TIMER_REGISTRAR_H
