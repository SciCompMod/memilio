.. _performance-monitoring-cpp:

Performance monitoring
======================

LIKWID
------

To measure the performance of the code, we use LIKWID. We recommend measuring the performance on an HPC system with LIKWID installed or refer to `<https://github.com/RRZE-HPC/likwid/wiki/Build>`_ on how to get LIKWID.
Run ``likwid-perfctr ./bin/my_example`` to measure the performance of the entire example. To measure the performance of a specific part, use the LIKWID marker API:

.. code-block:: cpp

    #include <likwid-marker.h>

    // ...

    LIKWID_MARKER_INIT;

    LIKWID_MARKER_START("my_marker");
    // code to be measured
    LIKWID_MARKER_STOP("my_marker");

    LIKWID_MARKER_CLOSE;

Set the CMake variable ``MEMILIO_USE_LIKWID=ON`` to enable LIKWID support and run ``likwid-perfctr -m ./bin/my_example``.
For more details see the LIKWID documentation, available `here <https://github.com/RRZE-HPC/likwid/wiki/likwid-perfctr>`_.


Performance timers
------------------

Here we present MEmilio's own timing framework. 


Timer usage
~~~~~~~~~~~

In this section we present how to use the AutoTimer class. This is the preferred way of using the timing framework, as
the class takes care of running the timer, managing its lifetime for later evaluation, and ensuring thread safety with
OpenMP. The AutoTimer class uses some other classes that are listed and explained in :ref:`Classes and their responsibilities`.

An AutoTimer starts when it is created, and stops when it is destroyed - which usually happens at the next closing
bracket :code:`}` or the next :code:`return`. This design, automating the starting and stopping of a timer, is
intentionally limiting, because it helps to avoid several issues or mistakes that can arise when manually running
timers. 

An example on how to use the timers can be found at `code example <https://github.com/SciCompMod/memilio/tree/main/cpp/examples>`__.


Timing in executables
^^^^^^^^^^^^^^^^^^^^^

To measure how long advancing a simulation without the setup takes, you can write 

.. code-block:: cpp

    #include "memilio/timer/auto_timer.h"

    int main() {
        Simulation sim; ... // setup

        {
            mio::timing::AutoTimer<"my simulation"> my_timer; // my_timer starts here
            sim.advance(t_max); // run the simulation
        } // my_timer stops here

        ... // evaluate results
    }

and will see a table printed at the end of your program next to the timer named "my simulation", that lists the time it
took to :code:`advance`.

You can add more timers like this, but make sure you use unique names, otherwise the same timer will be reused, and the
measured times will be added together. The name of the timer object itself (here :code:`my_timer`) is not important, as
long as the compiler does not complain about it.

Note that the automatic print can be disabled, and a manual print can be performed instead. There is also a print method
that gathers timers from all ranks in an MPI parallel context. Check out :code:`mio::timing::TimerRegistrar` for details
on these methods.


Timing in the library
^^^^^^^^^^^^^^^^^^^^^

Adding timers in the library is not much different to adding timers in main, but avoiding name collisions can be more
difficult. Hence, we use the optional second template argument of AutoTimer to specify its scope, as shown in the
following examples.

To measure the time a class member function takes, add a timer like this:

.. code-block:: cpp

    #include "memilio/timer/auto_timer.h"

    namespace foo

    class Bar {
        void baz() {
            mio::timing::AutoTimer<"baz", "foo::Bar"> timer;
            
            ... // rest of the function
        }
    };
    
    } // namespace foo

Or, when timing a free function:

.. code-block:: cpp

    #include "memilio/timer/auto_timer.h"

    namespace foo {

    void bar() {
        AutoTimer<"bar", "foo"> timer;
            
        ... // rest of the function
    }

    } // namespace foo

The first string given to AutoTimer is the timer's name, the second the scope it is in. They are used in combination
to identify the timer, similar to a map key, so they must be unique. This can be effectively guaranteed, if the name
matches the function and the scope contains all enclosing namespaces, like in the examples above.

If the containing function is used, a summary with timing results will be printed where both timers will show up as
:code:`foo::Bar::baz` and :code:`foo::bar`, respectively.


General recommendations
~~~~~~~~~~~~~~~~~~~~~~~

- **Do not time every detail.**
  While accurate, the timers are not very precise, so if you want to know how much time one or a few instructions take,
  use a profiler like (g)perf or likwid. Also, adding too many timers will clutter the timing results.

- **Only time computationally intensive code.**
  Similar to the last point, avoid timing small functions like setters and getters, and reserve using timers for the
  main compute loops or functions. While the timers add only a little overhead, it will become measurable when used too
  often.

- **Time entire functions.**
  Adding scopes for timing parts of main is fine, but you should avoid segmenting functions, either with scopes for
  AutoTimer or with manually run timers. The reason for this is related less to timers and more to code design, because
  if you can segment the function into multiple distinct parts, it is probably doing too many things, and should be
  separated into smaller functions. Also, adding scope (and thus indents) for AutoTimer does make code slightly harder
  to read.


The timing framework
~~~~~~~~~~~~~~~~~~~~

The main goals of this timing framework are to provide time measuring capabilities with minimal runtime overhead and
without having to plan around them. This means that accessing, starting and stopping a timer should be as fast as
possible, while the interfaces of the classes or functions that are to be timed should not change. Additionally, the
timer should work in parallel environments.

The solution to this is AutoTimer, whose usage was already shown above. There are, of course, some drawbacks. For
example, NamedTimer (the class used by AutoTimer) cannot be instantiated dynamically, as their name (and scope) have to
be known at compile time. This also means that adding a lot of timers will impact the time it takes to compile the code,
though a couple hundred timers should only take around an additional second.

.. _Classes and their responsibilities:
Classes and their responsibilities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this section, we describe the main components of the timing framework and how they interact. For more details on a
specific component, view its API documentation. 

- **BasicTimer**:
  The foundation of the timing framework. BasicTimer is a very simple class, that defines the methods start, stop,
  reset, and get_elapsed_time. These are used by all other classes in this framework. It uses a wall clock, so if compute
  resources are shared with other tasks, the timing results may be higher than expected. In debug builds, it will log
  errors whenever a member function was used incorrectly, e.g. when start was called twice.

- **TimerRegistration**:
  This simple struct is used to keep track of timers and some additional information, but does not manage their storage.
  It consists of two strings for name and scope, a reference to a BasicTimer, and a thread id. The thread id specifies
  which thread the timer is used in, which could differ from the thread it is created by.

- **Printer**:
  A pure virtual class defining a print method to evaluate and output timing results via a list of TimerRegistrations.
  Implemented by TablePrinter and ListPrinter.

- **TimerRegistrar**:
  Keeps track of timers via a list of TimerRegistrations, and holds a Printer that can be used to display all
  registered timers after the end of main. Timers can be registered by passing a TimerRegistration to its add_timer
  method. Uses a singleton pattern to provide global access to the same object, that is, the only way to obtain a
  TimerRegistrar object is by using its get_instance method, which returns a reference to a static object. Importantly,
  this class does not manage or own timer objects, and there is intentionally no methods that retrieve or delete
  TimerRegistrations.

- **NamedTimer**:
  Inherits from BasicTimer, with the main purpose of managing the lifetime, access, and registration of a timer.
  This is done using a singleton pattern, similar to TimerRegistrar, but the reference returned by get_instance is
  thread_local as well as static. The template parameters Name and Scope allow using more than one NamedTimer, since
  different template arguments define a different type. This effectively creates a global compile-time map, mapping a
  Name and Scope to a BasicTimer. Additionally, the NamedTimer registers itself automatically, and will only be
  destroyed after the TimerRegistrar.

- **AutoTimer**:
  Automates running an existing timer, by calling start in its constructor, and stop in its destructor. The timer used
  can be either specified via the Name and Scope template, fetching the corresponding NamedTimer internally, or by
  passing an lvalue reference to a BasicTimer.

Using NamedTimer and BasicTimer
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Preferably, you should use AutoTimer where possible, as its limiting design helps to avoid common errors, for example
with parallel regions. But, if you have to, you can use a NamedTimer directly without any extra work:

.. code-block:: cpp

    #include "memilio/timer/named_timer.h"

    namespace foo

    class Bar {
        void baz() {
            mio::timing::NamedTimer<"baz", "foo::Bar">::get_instance().start();
            
            ... // rest of the function

            mio::timing::NamedTimer<"baz", "foo::Bar">::get_instance().stop();
        }
    };
    
    } // namespace foo

This will behave exactly like the AutoTimer in the example above, while also allowing you to use the reset or
get_elapsed_time methods defined by BasicTimer.

Last but not least, you can also use a BasicTimer directly. This means that you will have to manually take care of
the timer object, threading and evaluation. If you add such a BasicTimer to the TimerRegistrar, you will probably need
to disable the final timer summary, and call print manually. Of course, you can also make your own list of registrations
and use a Printer directly.




Agent-based model benchmarks
----------------------------

There is a suite of benchmarks for the ABM that are used to check its performance. The suite contains setups of different sizes, to check that the model maintains its linear scaling. 

When you make any changes to the ABM or code used by it, run the benchmarks to check that its performance did not degrade. What exactly impacts the ABM's performance can be hard to tell (even parameter values may change its runtime), so it is best to run the bencharks on any change to the main library or the ABM specific code.

If you added a new feature (i.e., you didn't just fix a bug in an existing feature), make sure the feature is actually used by the benchmark. Add it to the benchmark if necessary, then run the benchmark to see if the cost for the new feature is acceptable and as expected.

Most new features will add some overhead, but this needs to be limited and in proportion to the added value of the feature so runtime doesn't grow out of control. Optional features that can be disabled should only incur minimal overhead. Always make sure there are no major performance regressions compared to the code in the current *main* branch.

Build the benchmarks by defining the CMake variable ``MEMILIO_BUILD_BENCHMARKS=ON`` in the build. Make sure to use a **Release** build to test performance.

.. code-block:: bash

    cmake .. -DMEMILIO_BUILD_BENCHMARKS=ON -DCMAKE_BUILD_TYPE=Release
    cmake --build .

Run the benchmark executable:

.. code-block:: bash

    ./build/bin/abm_benchmark

Each benchmark is run for a number of iterations and the average time is reported.

.. code-block:: text

    Benchmark                                 Time             CPU   Iterations
    ---------------------------------------------------------------------------
    abm_benchmark/abm_benchmark_50k        7583 ms         7583 ms            1
    abm_benchmark/abm_benchmark_100k      18216 ms        18214 ms            1
    abm_benchmark/abm_benchmark_200k      41492 ms        41489 ms            1

You may get a warning:

.. code-block:: text

    ***WARNING*** CPU scaling is enabled, the benchmark real time measurements may be noisy and will incur extra overhead.

If possible, disable CPU scaling to improve the consistency of results. See the Google Benchmark documentation here:
https://google.github.io/benchmark/reducing_variance.html

Also, try to reduce other system load during the benchmark run.

If it is not possible to disable frequency scaling, increase the runtime of the benchmark using the commands below. Constant CPU frequency is necessary to get the most reliable results and to measure small differences.

**REMINDER:** Don't forget to re-enable CPU scaling after you ran the benchmarks to save energy. Rebooting may restore the settings as well.

The benchmark executable has a number of command line arguments that customize execution. Use ``--help`` to list them all.

Two important options for consistency and stability:

- ``--benchmark_min_time=<T>``: Iterate each benchmark so that the total runtime is at least ``T`` seconds.  
  Default is 1 second, which may not be enough.  
  Try 60 seconds for better stability (you may need to experiment).

- ``--benchmark_repetitions=<N>``: Repeat each benchmark ``N`` times and report mean, median, and variance.  
  (Repetitions are **not** iterations; a benchmark can be repeated 10 times with 5 iterations each. Each repetition runs for at least the minimum time.)

``benchmark_repetitions`` is useful to check timing consistency, as it reports variance.  
However, it can be expensive because long-running benchmarks are repeated.  
``benchmark_min_time`` increases iterations only for fast-running benchmarks, which tend to be less stable.

**Suggested workflow:**

1. Use the benchmark to check the performance of your current changes:
  
  1. Run with 5–10 repetitions to check variance.
  2. Increase ``benchmark_min_time`` until variance is acceptable.
  3. Continue benchmarking with 1 repetition and the adjusted minimum time.

2. Repeat the benchmark *on the main branch* with the same ``benchmark_min_time`` and compare the results.
