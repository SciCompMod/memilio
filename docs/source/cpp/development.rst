Development
===========

.. _performance-monitoring-cpp:

Performance monitoring
----------------------

LIKWID
~~~~~~

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

Agent-based model benchmarks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
