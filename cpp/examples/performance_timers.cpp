#include "memilio/timer/definitions.h"
#include "memilio/timer/table_printer.h"
#include "memilio/timer/timer_registrar.h"
#include "memilio/utils/mioomp.h"
#include "memilio/timer/timers.h"

#include <thread>

/// @brief Workload for the timers to report on.
void load()
{
    // Two very simple work loads to test the timers. Choose one, an comment out the other.

    // Option a: increment a counter. Extremely cheap load, so for large N, the timing overhead will matter.
    // Uses volatile in an attempt to avoid compiler optimisations.

    // volatile static int ctr = 0;
    // ++ctr;

    // Option b: let the current thread briefly wait. The timing overhead should be negligible.

    std::this_thread::sleep_for(std::chrono::milliseconds(3));
}

int main()
{
    // Start immediately timing the main function. An AutoTimer starts timing upon its creation, and ends timing when
    // it is destroyed. This usually happens when a function returns, or a scope indicated by curly braces "{}" ends.
    // The name of the AutoTimer does not matter, but the name given to its constructor is used to identify which
    // NamedTimer to use.
    mio::AutoTimer timer_ms(CONST_LITERAL("main"));

    // Set a Printer to output the timers.
    // These only have to be set if you want to customize the output. By default, TimerRegistrar prints all timers at
    // the end of the programm. This can be disabled by calling `TimerRegistrar::disable_final_timer_summary()`.
    auto printer = std::make_unique<mio::timing::TablePrinter>();
    printer->set_time_format("{:e}");
    mio::timing::TimerRegistrar::get_instance().set_printer(std::move(printer));

    // To manually print all timers, use `TimerRegistrar::print_timers()`, but make sure that no timers are running.
    // The "main" timer in this example would make that difficult, but you can simply add another scope.

    const int N = 1000; // Number of iterations. Be wary of increasing this number when using the sleep_for load!

    std::cout << "Num threads: " << mio::get_omp_max_threads() << "\n";

    // Open a new scope to time computations.
    {
        mio::AutoTimer timer_cs(CONST_LITERAL("computation scope"));

        // First compute loop.
        {
            mio::AutoTimer timer_out(CONST_LITERAL("first compute loop (without inner timer)"));

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                load();
            }
        }
        // Second compute loop, with additional timers during computation. While the timing overhead is minimal, it is
        // measurable, so it should be generally avoided to time inside main compute loops.
        {
            mio::AutoTimer timer_out(CONST_LITERAL("second compute loop (with inner timer)"));

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                mio::AutoTimer timer_in(CONST_LITERAL("  compute loop inner timer"));
                load();
            }
        }
    }

    return 0;
}
