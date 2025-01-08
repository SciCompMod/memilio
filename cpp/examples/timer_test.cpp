#include "memilio/utils/mioomp.h"

#include <chrono>
#include <ctime>
#include <iostream>
#include <list>
#include <map>
#include <mutex>
#include <string>
#include <thread>
#include <utility>

// #define DEBUG(cout_args) std::cerr << cout_args << std::endl << std::flush;
#define DEBUG(cout_args)

//////////////////////
// src: https://stackoverflow.com/questions/15858141/conveniently-declaring-compile-time-strings-in-c/15863804#15863804
//////////////////////

// helper function
constexpr unsigned c_strlen(char const* str, unsigned count = 0)
{
    return ('\0' == str[0]) ? count : c_strlen(str + 1, count + 1);
}

// destination "template string" type
template <char... chars>
struct exploded_string {
    inline static constexpr char data[] = {chars...};
    static void print()
    {
        char const str[] = {chars...};
        std::cout.write(str, sizeof(str));
    }
};

// struct to explode a `char const*` to an `exploded_string` type
template <typename StrProvider, unsigned len, char... chars>
struct explode_impl {
    using result = typename explode_impl<StrProvider, len - 1, StrProvider::str()[len - 1], chars...>::result;
};

// recursion end
template <typename StrProvider, char... chars>
struct explode_impl<StrProvider, 0, chars...> {
    using result = exploded_string<chars...>;
};

// syntactical sugar
template <typename StrProvider>
using explode = typename explode_impl<StrProvider, c_strlen(StrProvider::str())>::result;

//////////////////////
//////////////////////
//////////////////////

namespace timing
{

using ClockType    = std::chrono::steady_clock;
using TimeType     = ClockType::time_point;
using DurationType = ClockType::duration;

// TODO: add thread id?
using TimerRegistration = std::pair<const char*, const DurationType*>;

inline double printable_time(DurationType time)
{
    return std::chrono::duration_cast<std::chrono::duration<double>>(time).count();
}

template <const char Name[]>
struct my_str_provider {
    constexpr static char const* str()
    {
        return Name;
    }
};

template <char... Chars>
class Timer2;

template <char... Chars>
Timer2<Chars...> timer_type(exploded_string<Chars...>);

template <const char Name[]>
using timer_from_name = decltype(timer_type(explode<my_str_provider<Name>>{}));

// holds raw pointers to names and durations of the thread_local Timer instances
class TimerRegistrar
{
public:
    void add_timer(const TimerRegistration& registration)
    {
        m_registration_lock.lock();
        m_registry.push_back(registration);
        m_registration_lock.unlock();
    }

    // preregistration does not seem to affect timing basically at all. it removes the registration from timers started
    // before a preregistered timer, but a timer does not measure its own registration
    //
    // TODO: does this need to be in a omp parallel to correctly preregister?
    template <const char*... Names>
    void preregister_timers()
    {
        (timer_from_name<Names>::get_instance(), ...);
    }

    // questionable whether this should be allowed to be used over print_on_exit.
    void print_timers(std::ostream& out = std::cout)
    {
        PRAGMA_OMP(single)
        {
            const auto indent = "  ";
            out << "All Timers: " << m_registry.size() << "\n";
            for (const auto& entry : m_registry) {
                out << indent << entry.first << ": " << std::scientific << printable_time(*entry.second) << "\n";
            }
            // dedupe list entries from parallel execution
            std::map<std::string, DurationType> deduper;
            for (const auto& entry : m_registry) {
                deduper[entry.first] = DurationType::zero(); // is this neccessary?
            }
            for (const auto& entry : m_registry) {
                deduper[entry.first] += *entry.second;
            }
            out << "Unique Timers (accumulated): " << deduper.size() << "\n";
            for (const auto& entry : deduper) {
                out << indent << entry.first << ": " << std::scientific << printable_time(entry.second) << "\n";
            }
        }
    }

    void print_on_exit() // requires dtor (i.e. get_instance scheme)
    {
        m_print_on_death = true;
    }

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

private:
    TimerRegistrar()
    {
        DEBUG("### Registrar at your service")
    }

    ~TimerRegistrar()
    {
        if (m_print_on_death) {
            PRAGMA_OMP(single)
            {
                std::cout << "Final Timer Summary\n";
                print_timers();
            }
            DEBUG("### Registrar retired")
        }
    }

    std::list<TimerRegistration> m_registry;
    std::mutex m_registration_lock;
    bool m_print_on_death;
};

// since it is thread_local and typed by its name, no collisions can occur between different timers or threads
//
// though one can still double start or stop the timer, creating bogus results. this can be fixed by making everything
// protected and using e.g. ScopedTimer with inheritance
//
// BUG: timers can have same name, but different type, since their type differs by pointer...
// alternatives: Timer2, enums, C++20
template <const char Name[]>
class Timer
{
public:
    static Timer& get_instance()
    {
        thread_local static Timer t;
        return t;
    }

    void start()
    {
        m_start_time = ClockType::now();
    }

    void stop()
    {
        const TimeType stop_time = ClockType::now();
        m_elapsed_time += stop_time - m_start_time;
    }

private:
    Timer()
        : m_start_time()
        , m_elapsed_time(DurationType::zero())
    {
        register_self();
        DEBUG("--- Timer registration : " << Name)
    }

    ~Timer()
    {
        DEBUG("--- Timer destruction  : " << Name)
    }

    void register_self() // inline into ctor?
    {
        TimerRegistrar::get_instance().add_timer({Name, &m_elapsed_time});
    }

    TimeType m_start_time;
    DurationType m_elapsed_time;
};

// second attempt, templates dont look nice at all, but are indep. of the char[] used, only using its content
template <char... Chars>
class Timer2
{
public:
    static constexpr auto Name = exploded_string<Chars...>::data;

    static Timer2& get_instance()
    {
        thread_local static Timer2 t;
        return t;
    }

    void start()
    {
        m_start_time = ClockType::now();
    }

    void stop()
    {
        const TimeType stop_time = ClockType::now();
        m_elapsed_time += stop_time - m_start_time;
        DEBUG("    Elapsed: " << printable_time(m_elapsed_time))
    }

private:
    Timer2()
        : m_start_time()
        , m_elapsed_time(DurationType::zero())
    {
        register_self();
        DEBUG("--- Timer registration : " << Name)
    }

    ~Timer2()
    {
        DEBUG("--- Timer destruction  : " << Name)
    }

    void register_self() // inline into ctor?
    {
        TimerRegistrar::get_instance().add_timer({Name, &m_elapsed_time});
    }

    TimeType m_start_time;
    DurationType m_elapsed_time;
};

// specialization allowing disabling timers by changing their name. good method?
// alternative : constexpr ifs? keeps stuff in place, but bloates Timer2
template <>
class Timer2<'@', 'd', 'i', 's', 'a', 'b', 'l', 'e'>
{
public:
    static constexpr auto Name = "@disable";

    static Timer2& get_instance()
    {
        thread_local static Timer2 t;
        return t;
    }

    void start()
    {
    }

    void stop()
    {
    }
};

// starts on creating, stops at end of scope. can be called repeatedly
//
// potential misuse by starting this too late in a scope
template <const char Name[]>
class ScopedTimer
{
public:
    ScopedTimer()
    {
        timer_from_name<Name>::get_instance().start();
        DEBUG("... ScopedTimer start: " << Name)
    }
    ~ScopedTimer()
    {
        DEBUG("... ScopedTimer stop:  " << Name)
        timer_from_name<Name>::get_instance().stop();
    }
};

} // namespace timing

int main()
{
    DEBUG("*** enter main")

    static constexpr char timer_name_hw[]                           = "@disable";
    static constexpr char timer_name_ms[]                           = "main scope";
    static constexpr char timer_name_loop_inner[]                   = "Loop Inner";
    static constexpr char timer_name_loop_inner_prereg[]            = "Loop Inner Prereg";
    static constexpr char timer_name_loop_outer[]                   = "Loop Outer no Inner";
    static constexpr char timer_name_loop_outer_with_inner[]        = "Loop Outer with Inner";
    static constexpr char timer_name_loop_outer_with_inner_prereg[] = "Loop Outer with Inner Prereg";
    timing::ScopedTimer<timer_name_ms> timer_ms;

    timing::TimerRegistrar::get_instance().preregister_timers<timer_name_loop_inner_prereg>();
    timing::TimerRegistrar::get_instance().print_on_exit();

    volatile int ctr = 0;
    int N            = 10000;
    (void)ctr;
    const auto load = [&] {
        ++ctr;
        // std::this_thread::sleep_for(std::chrono::milliseconds(3));
    };

    DEBUG("*** main post setup")
    std::cout << "Num threads: " << mio::get_num_threads() << "\n";

    {
        timing::ScopedTimer<timer_name_hw> timer_hw;

        {
            timing::ScopedTimer<timer_name_loop_outer_with_inner> timer_loop_outer_with_inner;

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                timing::ScopedTimer<timer_name_loop_inner> timer_loop_inner;
                load();
            }
        }
        {
            timing::ScopedTimer<timer_name_loop_outer_with_inner_prereg> timer_loop_outer_with_inner;

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                timing::ScopedTimer<timer_name_loop_inner_prereg> timer_loop_inner;
                load();
            }
        }

        {
            timing::ScopedTimer<timer_name_loop_outer> timer_loop_outer;

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                load();
            }
        }
    }

    timing::TimerRegistrar::get_instance().print_timers();

    DEBUG("*** exit main")
}
