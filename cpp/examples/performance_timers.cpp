#include "memilio/timer/list_printer.h"
#include "memilio/timer/table_printer.h"
#include "memilio/timer/timer_registrar.h"
#include "memilio/utils/mioomp.h"
#include "memilio/timer/timers.h"

#include <functional>
#include <thread>

int main()
{
    mio::ScopedTimer timer_hw(CONST_LITERAL("Hello World!"));

    // Set a Printer to output the timers.
    // Both printer and format used here are the default value, we set them only as a demonstration how they could be
    // changed. If the defaults are
    auto tf = std::make_unique<mio::timing::TablePrinter>();
    tf->set_time_format("{:e}");
    mio::timing::TimerRegistrar::get_instance().set_printer(std::move(tf));

    volatile int ctr = 0;
    int N            = 1000;
    (void)ctr;
    const auto load = [&] {
        // ++ctr;
        std::this_thread::sleep_for(std::chrono::milliseconds(3));
    };

    std::cout << "Num threads: " << mio::get_omp_max_threads() << "\n";

    {
        mio::ScopedTimer timer_ms(CONST_LITERAL("main"));

        {
            mio::ScopedTimer timer_out(CONST_LITERAL("outside loop (with inner timer)"));

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                mio::ScopedTimer timer_in(CONST_LITERAL("inside loop"));
                load();
            }
        }
        {
            mio::ScopedTimer timer_out(CONST_LITERAL("outside loop"));

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                load();
            }
        }

        {
            mio::timing::get_timer(CONST_LITERAL("outside loop (manual timer)")).start();

            PRAGMA_OMP(parallel for)
            for (int i = 0; i < N; i++) {
                load();
            }

            mio::timing::get_timer(CONST_LITERAL("outside loop (manual timer)")).stop();
        }
    }

    return 0;
}
