/* 
* Copyright (C) 2020-2026 MEmilio
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
#include "memilio/timer/auto_timer.h"
#include "memilio/timer/table_printer.h"

#include <thread> // This is only used for the example load function.

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
    // Specify the namespace of AutoTimer, so we don't have to repeat it. Do not do this with entire namespaces.
    using mio::timing::AutoTimer;

    // Start immediately timing the main function. An AutoTimer starts timing upon its creation, and ends timing when
    // it is destroyed. This usually happens when a function returns, or a scope indicated by {curly braces} ends.
    // The name of the AutoTimer object (here timer_ms) does not matter, but the name given to it as a template
    // (i.e. "main") does, as it is used to identify the NamedTimer that is used internally, it works like a map key.
    AutoTimer<"main"> timer_ms;

    // Set a Printer to output the timers.
    // These only have to be set if you want to customize the output. By default, TimerRegistrar prints all timers at
    // the end of the programm. This can be disabled by calling `TimerRegistrar::disable_final_timer_summary()`.
    auto printer = std::make_unique<mio::timing::TablePrinter>();
    printer->set_time_format("{:e}");
    mio::timing::TimerRegistrar::get_instance().set_printer(std::move(printer));

    // To manually print all timers, use `TimerRegistrar::print_timers()`, but make sure that no timers are running.
    // The "main" timer in this example would make that difficult, but you can simply add another scope around it,
    // similar to the "compute loops" timer below.

    const int N = 1000; // Number of iterations. Be wary of increasing this number when using the sleep_for load!

    std::cout << "Num threads: " << mio::omp::get_max_threads() << "\n";

    // Open a new scope to time computations.
    {
        // This is another a good example of how to use AutoTimer, now inside an unnamed scope.
        AutoTimer<"compute loops"> timer_cs;

        // First compute loop.
        {
            // With the second string we define the scope of the timer. Here it is kind of misused, the scope should
            // usually name the surrounding class, struct, or function. This has two reasons: First, it makes name
            // collisions less likely. Second, the scope can be used to analyze the results. Hence, "first loop" would
            // be way too granular, the appropriate scope in this case would be "main", if any.
            AutoTimer<"compute", "first loop"> timer_out;

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                load();
            }
        }
        // Second compute loop, with an additional timer during computation. While the timing overhead is minimal, it is
        // still measurable, so it should be generally avoided to time inside main compute loops.
        {
            // Again, we use a not quite optimal name for the scope, but it is shared with the inner timer. In this
            // example, it would be enough to include "second loop" in the name. The scope should be used for
            // differentiating between, e.g., different models or simulations.
            AutoTimer<"compute (with inner timer)", "second loop"> timer_out;

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                // Note that when using multithreading, this timer will measure a longer (aggregated) time than its
                // outer timer.
                AutoTimer<"inner timer", "second loop"> timer_in;
                load();
            }
        }
    }

    return 0;
}
