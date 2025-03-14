#include "memilio/utils/mioomp.h"
#include "memilio/utils/stl_util.h"

#include <chrono>
#include <ctime>
#include <iostream>
#include <list>
#include <map>
#include <mutex>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

// clang-format off
#define CONST_LITERAL(literal) []() constexpr { return literal; }
// clang-format on

namespace mio
{

namespace details
{

constexpr std::size_t constexpr_strlen(const char* string_literal)
{
    // returns length of string, does not count terminating 0 char
    return string_literal[0] == '\0' ? 0 : constexpr_strlen(string_literal + 1) + 1;
}

} // namespace details

namespace timing
{

#ifdef MEMILIO_ENABLE_OPENMP

using TimeType     = decltype(omp_get_wtime());
using DurationType = decltype(omp_get_wtime());

inline double time_in_seconds(DurationType time)
{
    return time;
}

inline TimeType get_time_now()
{
    return omp_get_wtime();
}

#else

using TimeType     = std::chrono::steady_clock::time_point;
using DurationType = std::chrono::steady_clock::duration;

inline double time_in_seconds(DurationType time)
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(time).count();
}

inline TimeType get_time_now()
{
    return std::chrono::steady_clock::now();
}

#endif

class BasicTimer
{
public:
    void start()
    {
        m_start_time = get_time_now();
    }

    void stop()
    {
        const TimeType stop_time = get_time_now();
        m_elapsed_time += stop_time - m_start_time;
    }

    DurationType get_elapsed_time() const
    {
        return m_elapsed_time;
    }

private:
    TimeType m_start_time;
    DurationType m_elapsed_time;
};

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

// second attempt, templates dont look nice at all, but are indep. of the char[] used, only using its content
template <class Tag, char... Chars>
class Timer : public BasicTimer
{
public:
    static constexpr char Name[] = {Chars..., '\0'};

    static Timer& get_instance()
    {
        // create and return a timer
        // the constructor is only called the first time this function is used (for each thread)
        thread_local static Timer t;
        return t;
    }

private:
    Timer()
    {
        TimerRegistrar::get_instance().add_timer({Name, *this, mio::get_omp_thread_id()});
    }
};

namespace details
{

template <class Tag, typename Identifier, std::size_t... I>
constexpr auto& get_timer_by_name(Identifier timer_name, std::index_sequence<I...>)
{
    (void)timer_name; // if the name is "", i.e. empty, the Identifier shows up as unused
    return Timer<Tag, timer_name()[I]...>::get_instance();
}

} // namespace details

template <class Tag = void, class Identifier>
constexpr auto& get_timer(Identifier timer_name)
{
    static_assert(std::is_same_v<decltype(std::declval<Identifier>()()), const char*>, "message");
    return details::get_timer_by_name<Tag>(timer_name,
                                           std::make_index_sequence<mio::details::constexpr_strlen(timer_name())>{});
}

template <class Tag = void>
class ScopedTimer
{
public:
    template <class Identifier>
    ScopedTimer(Identifier timer_name)
        : m_timer(timing::get_timer<Tag>(timer_name))
    {
        m_timer.start();
    }

    ScopedTimer(timing::BasicTimer& timer)
        : m_timer(timer)
    {
        m_timer.start();
    }

    ~ScopedTimer()
    {
        m_timer.stop();
    }

private:
    timing::BasicTimer& m_timer;
};

} // namespace timing

using timing::ScopedTimer;

} // namespace mio

int main()
{
    mio::ScopedTimer timer_hw(CONST_LITERAL("Hello World!"));

    // timing::TimerRegistrar::get_instance().preregister_timers<timer_name_loop_inner_prereg>();
    mio::timing::TimerRegistrar::get_instance().print_on_exit();

    volatile int ctr = 0;
    int N            = 1000;
    (void)ctr;
    const auto load = [&] {
        // ++ctr;
        std::this_thread::sleep_for(std::chrono::milliseconds(3));
    };

    std::cout << "Num threads: " << mio::get_omp_num_threads() << "\n";

    {
        mio::ScopedTimer timer_ms(CONST_LITERAL("Main"));

        {
            mio::ScopedTimer timer_out(CONST_LITERAL("outside loop with inner timer"));

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                mio::ScopedTimer timer_in(CONST_LITERAL("inside loop"));
                load();
            }
        }
        {
            mio::ScopedTimer timer_out(CONST_LITERAL(""));

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                load();
            }
        }

        {
            auto& timer_out = mio::timing::get_timer(CONST_LITERAL("outside loop manual timer"));
            timer_out.start();

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                load();
            }

            timer_out.stop();
        }
    }

    mio::timing::TimerRegistrar::get_instance().print_timers();
}
