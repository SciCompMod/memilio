#include "memilio/utils/mioomp.h"
#include <cassert>
#include <chrono>
#include <cstddef>
#include <ctime>
#include <iostream>
#include <list>
#include <map>
#include <mutex>
#include <numeric>
#include <string>
#include <string_view>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

// #define DEBUG(cout_args) std::cerr << cout_args << std::endl << std::flush;
#define DEBUG(cout_args)

namespace timing
{

#ifdef MEMILIO_ENABLE_OPENMP

using TimeType     = omp_wtime_t;
using DurationType = omp_wtime_t;

inline double time_in_seconds(DurationType time)
{
    return time;
}

inline TimeType time_now()
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

inline TimeType time_now()
{
    return std::chrono::steady_clock::now();
}

#endif

// struct TimerRegistration {
//     const std::string name;
//     const std::vector<std::string> scope;
//     DurationType elapsed_time;
//     const int thread_id;
// };

class Timer
{
public:
    Timer(const std::string& name, const std::vector<std::string>& scope)
        : m_name(name)
        , m_scope(scope)
        , m_start_time()
        , m_elapsed_time()
    {
    }

    void start()
    {
        m_start_time = time_now();
    }

    void stop()
    {
        const TimeType stop_time = time_now();
        m_elapsed_time += stop_time - m_start_time;
    }

    // const std::string_view name() const
    // {
    //     return m_name;
    // }

    DurationType get_elapsed_time() const
    {
        return m_elapsed_time;
    }

private:
    std::string m_name;
    std::vector<std::string> m_scope;
    TimeType m_start_time;
    DurationType m_elapsed_time;
};

// holds raw pointers to names and durations of the thread_local Timer instances
class TimerRegistrar
{
public:
    // void add_timer(const std::string& name, const std::vector<std::string>& scope)
    // {
    //     m_registration_lock.lock();
    //     m_registry[mio::get_omp_thread_id()].emplace(name, scope);
    //     m_registration_lock.unlock();
    // }

    Timer& get_timer(const std::string& name, const std::vector<std::string>& scope)
    {
        const int thread_id = mio::get_omp_thread_id();
        auto timer_itr      = m_registry[thread_id].find(name);
        if (timer_itr == m_registry[thread_id].end()) {
            m_registration_lock.lock();
            const auto emplace_result = m_registry[thread_id].emplace(name, Timer{name, scope});
            assert(emplace_result.second && "Bug: find says the timer does not exist, emplace says it does.");
            timer_itr = emplace_result.first;
            m_registration_lock.unlock();
        }
        return timer_itr->second;
    };

    // questionable whether this should be allowed to be used over print_on_exit.
    void print_timers(std::ostream& out = std::cout)
    {
        PRAGMA_OMP(single)
        {
            const auto indent = "  ";
            out << "All Timers: " << std::accumulate(m_registry.begin(), m_registry.end(), 0, [](auto&& s, auto&& reg) {
                return s + reg.size();
            }) << "\n";
            for (unsigned thread_id = 0; thread_id < m_registry.size(); thread_id++) {
                for (const auto& [name, timer] : m_registry[thread_id]) {
                    out << indent << name << ": (" << thread_id << ")" << std::scientific
                        << time_in_seconds(timer.get_elapsed_time()) << "\n";
                }
            }
            // dedupe list entries from parallel execution
            std::map<std::string, DurationType> deduper;
            for (unsigned thread_id = 0; thread_id < m_registry.size(); thread_id++) {
                for (const auto& [name, timer] : m_registry[thread_id]) {
                    deduper[name] += timer.get_elapsed_time();
                }
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
        : m_registry(mio::get_omp_num_threads())
        , m_registration_lock()
        , m_print_on_death(false)
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

    std::vector<std::unordered_map<std::string, Timer>> m_registry;
    std::mutex m_registration_lock;
    bool m_print_on_death;
};

// template <const char Name[]>
class ScopedTimer
{
public:
    ScopedTimer(const std::string& name, const std::vector<std::string>& scope = {})
        : m_timer(TimerRegistrar::get_instance().get_timer(name, scope))
    {
        m_timer.start();
        DEBUG("... ScopedTimer start: " << Name)
    }
    ~ScopedTimer()
    {
        DEBUG("... ScopedTimer stop:  " << Name)
        m_timer.stop();
    }

private:
    Timer& m_timer;
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
    timing::ScopedTimer timer_ms(timer_name_ms);

    // timing::TimerRegistrar::get_instance().preregister_timers<timer_name_loop_inner_prereg>();
    timing::TimerRegistrar::get_instance().print_on_exit();

    volatile int ctr = 0;
    int N            = 1000;
    (void)ctr;
    const auto load = [&] {
        // ++ctr;
        std::this_thread::sleep_for(std::chrono::milliseconds(3));
    };

    DEBUG("*** main post setup")
    std::cout << "Num threads: " << mio::get_omp_num_threads() << "\n";

    {
        timing::ScopedTimer timer_hw(timer_name_hw);

        {
            timing::ScopedTimer timer_loop_outer_with_inner(timer_name_loop_outer_with_inner);

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                timing::ScopedTimer timer_loop_inner(timer_name_loop_inner);
                load();
            }
        }
        {
            timing::ScopedTimer timer_loop_outer_with_inner(timer_name_loop_outer_with_inner_prereg);

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                timing::ScopedTimer timer_loop_inner(timer_name_loop_inner_prereg);
                load();
            }
        }

        {
            timing::ScopedTimer timer_loop_outer(timer_name_loop_outer);

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                load();
            }
        }
    }

    // timing::TimerRegistrar::get_instance().print_timers();

    DEBUG("*** exit main")
}
