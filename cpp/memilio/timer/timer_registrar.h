#ifndef MIO_TIMER_TIMER_REGISTRAR_H
#define MIO_TIMER_TIMER_REGISTRAR_H

#include "memilio/timer/registration.h"
#include "memilio/timer/table_printer.h"

#include <list>
#include <memory>
#include <mutex>
#include <ostream>
#include <string>
#include <utility>

namespace mio
{
namespace timing
{

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

    const auto& get_register() const
    {
        return m_register;
    }

    void print_timers(std::ostream& out = std::cout) const
    {
        PRAGMA_OMP(single)
        {
            get_instance().m_printer->print(m_register, out);
        }
    }

    void disable_final_timer_summary() const
    {
        get_instance().m_print_on_death = false;
    }

    void set_printer(std::unique_ptr<Printer>&& printer)
    {
        m_printer.swap(printer);
        // old value of m_printer (now stored in printer) is deleted at end of scope
    }

private:
    TimerRegistrar() = default;

    ~TimerRegistrar()
    {
        if (m_print_on_death) {
            PRAGMA_OMP(single)
            {
                std::cout << "Final Timer Summary\n";
                m_printer->print(m_register, std::cout);
            }
        }
    }

    std::unique_ptr<Printer> m_printer = std::make_unique<TablePrinter>();
    bool m_print_on_death              = true;
    std::list<TimerRegistration> m_register;
    std::mutex m_registration_lock;
};

} // namespace timing
} // namespace mio

#endif // MIO_TIMER_TIMER_REGISTRAR_H
