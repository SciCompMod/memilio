# Profiling

This document explains how to create performance profiles, i.e., measurements of where in the code time is spent, usually aggregated by function. There is a variety of open-source or proprietary tools available. Generally, there are three types of profilers:
- instrumenting: modify the compiled code to include code for measurement. Advantage: exact function call counts and times, no approximations. Disadvantage: instrumentation adds runtime overhead that can distort the result or disable compiler optimizations (e.g. functions aren't inlined), so short functions usually must be excluded.
- sampling: checking the CPU stack register in regular intervals to see which code is currently being executed. Advantage: usually low measurement runtime overhead. Disadvantage: stochastic results, imprecise especially for short functions, longer runtimes improve averages.
- emulating: running the program in a virtual environment. Advantage: exact and able to measure even small functions. Disadvantage: program can run significantly longer, no real world timing results but instruction cycle counts are generally directly translatable.

These tools have been tried for memilio and are described below:
- gperftools: sampling profiler that gives reliable results with little setup
- Score-P: large suite that supports both sampling and instrumenting. Special support for parallel code (OpenMP and MPI), including tracing to find load imbalances in the parallelization. 
- Valgrind: suite that includes the Callgrind emulating profiler; detailed descripton will be added soon, until then see the [official documentation](https://valgrind.org/info/tools.html).

How to use these tools, the general optimization loop:
1. Run a benchmark to get a baseline runtime and identify the need for performance optimization, e.g., a change introduced a performance regression.
2. Create a profile to see where time is lost; maybe compare with a profile before your changes.
3. Optimize the code.
4. Test the functionality.
5. Confirm the optimization with a benchmark and profile.
6. Repeat if necessary.

## gperftools

gperftools (formerly google performance tools) is a suite of performance tools. Here the focus is on the [profiler](https://gperftools.github.io/gperftools/cpuprofile.html). It's a small sampling profiler that is quick to setup and use on most systems. It's integrated into our build system for even greater convenience.

Basic steps:
1. Install gperftools. It's available in many package managers (apt packages `google-perftools` and `libgoogle-perftools-dev`, spack or homebrew package `gperftools`). E.g., on the DLR-SC hpda cluster (insert the compiler version you are using): `module load spack-user; module load PrgEnv/gcc13-openmpi; spack install gperftools%gcc@13.3.0; module load gperftools`.
2. configure memilio with cmake variable `MEMILIO_ENABLE_PROFILING=ON`. Note that compiling with profiling enabled does not incur any runtime overhead unless profling is also enabled at runtime (see step 4), so developers can just enable it always for convenience.
3. compile memilio.
4. run the program with environment variable `CPUPROFILE=profile.out` set.
5. generate a human-readable annotated call graph: `pprof-symbolize --pdf <exe> profile.out > profile.pdf`. Check the documentation for other output formats.

This generates a profile of the whole program. For each function, the graph contains the number of times the function was sampled and how much time is spent in the function as a percentage of the total runtime. Both sample count and percentage are displayed as exclusive numbers (i.e. the function without child functions) and inclusive numbers (the function including its child functions). 

Probably not all parts of the program are interesting to profile. It's possible to filter the profile, e.g., to look at specific functions, using the `--focus` flag of `pprof-symbolize` or by starting the profiling manually using the `MEMILIO_PROFILER_START(<filename>)` and `MEMILIO_PROFILER_STOP` macros in `memilio/utils/profiler.h`. E.g., when profiling the `abm_simulation` program, it may be a good idea to exclude pre- and postprocessing by starting profiling right before the call to `Simulation::advance`.

## Score-P

Score-P is a tool that can be used to measure the perfomance of parallel code and supports profiling and tracing. For information on how to get it, please have a look at the [Score-P quickstart](https://scorepci.pages.jsc.fz-juelich.de/scorep-pipelines/docs/scorep-6.0/html/quickstart.html). It can be installed with the spack manager (and probably many other package managers) or loaded as a module on hpc systems.
For the visualization of the performance data you need Cube (information on installation, see [here](https://apps.fz-juelich.de/scalasca/releases/cube/4.3/docs/CubeInstall.pdf)). It can load both Score-P profiles and Scalasca traces and you can find details about the loading and usage below or in the [User Guide](https://www.vi-hps.org/cms/upload/packages/cube/CubeGuide.pdf).

This file explains how to use Score-P with GCC. For other compilers this might look different or not work correctly.
Once you have Score-P you can build the abm with: 
```
export SCOREP_WRAPPER=off #disable scorep during cmake configure
cmake .. -DMEMILIO_ENABLE_OPENMP=ON -DCMAKE_CXX_COMPILER=${Path to scorep}/scorep-mpic++ -DCMAKE_C_COMPILER=${Path to scorep}/scorep-mpicc -DCMAKE_BUILD_TYPE=RelWithDebInfo -DMEMILIO_ENABLE_WARNINGS_AS_ERRORS=OFF
export SCOREP_WRAPPER=on #reenable scorep
export SCOREP_WRAPPER_INSTRUMENTER_FLAGS="--pomp --thread=omp --mpp=none --instrument-filter=${Path to filter}/scorep-filter-abm"
cmake --build . -j --target abm_simulation
```

It is useful to use ```--target``` to set the target to build, because Score-P has difficulties with some of the other targets.
The [filter file](./scorep-filter-abm) is designed to include only the interesting information about the abm. It filters (i.e. measures) only the basic framework of the abm, that is only the simulation part up to the loop over time. The measurement of small and frequently called functions increases the overhead of the profiling. This mainly concerns the runtime and makes the profile unreliable, since the time measured by the profiler is not the time that the functions would take in normal operation. To prevent this, those functions are neglected in the current filter and the profile may give a hint on which part of the abm you should investigate further. To measure additional functions you add ```INCLUDE *Class::myfunction*``` to the [filter](./scorep-filter-abm). 

To avoid the build errors about "include style" and "gcc extension" set ```-DMEMILIO_ENABLE_WARNINGS_AS_ERRORS=OFF``` during cmake configuration.

### Profile measurements

To run the basic analysis with profiling:
```
export OMP_NUM_THREADS=8 #or whatever
./bin/abm_simulation 1 <some result dir>
```

This gives a summary of how much time was spent in each function. There will be a directory named scorep_... with a .cubex file, that can be opened in the Cube visualization tool with ```cube scorep-folder/profile.cubex```. This can give a first indication for bottlenecks, then you don't need Scalasca.
Also, you can use ```scorep-score -r scorep-folder/profile.cubex``` to get a first impression. The output for a profile without any filter looks like this (shortened to a few functions):

    Estimated aggregate size of event trace:                   267GB
    Estimated requirements for largest trace buffer (max_buf): 267GB
    Estimated memory requirements (SCOREP_TOTAL_MEMORY):       267GB
    (warning: The memory requirements cannot be satisfied by Score-P to avoid
    intermediate flushes when tracing. Set SCOREP_TOTAL_MEMORY=4G to get the
    maximum supported memory or reduce requirements using USR regions filters.)

    flt     type      max_buf[B]         visits time[s] time[%] time/visit[us]  region
         ALL 286,150,801,259 11,005,733,830 2854.35   100.0           0.26  ALL
         USR 286,146,847,402 11,005,647,977 1626.05    57.0           0.15  USR
         OMP       3,931,924         85,010  152.08     5.3        1788.98  OMP
         COM          21,892            842 1076.21    37.7     1278156.85  COM
      SCOREP              41              1    0.01     0.0        5326.90  SCOREP

         USR 252,154,970,496  9,698,268,096 1292.41    45.3           0.13  mio::abm::Location&?mio::abm::Person::get_location()
        USR   7,055,253,036    271,355,886   40.67     1.4           0.15  void?Eigen::internal::ignore_unused_variable(const?T&)??with?T?=?long?int?
        USR   3,551,985,216    136,614,816   17.11     0.6           0.13  mio::abm::InfectionState?mio::abm::Person::get_infection_state(mio::abm::TimePoint)?const
        USR   2,106,780,026     81,030,001   23.04     0.8           0.28  size_t?mio::flatten_index(const?MultiIndex&,?const?MultiIndex&)??with?MultiIndex?=?mio::Index<mio::abm::VirusVariant,?mio::AgeGroup>?
    ...

When changing the abm one should check how much memory is used by Score-P and that it does not exceed 4 GB which is the maximum value for SCOREP_TOTAL_MEMORY. Also when tracing SCOREP_TOTAL_MEMORY should be adjusted to the value proposed by scorep-score. 

As you can see, the profile contains a lot of functions that are not of interest, for example small functions like ```get_location()``` or internal functions like ```ignore_unused_variables()```. Those functions are filtered when running with the proposed [filter](./scorep-filter-abm):

    Estimated aggregate size of event trace:                   369kB
    Estimated requirements for largest trace buffer (max_buf): 369kB
    Estimated memory requirements (SCOREP_TOTAL_MEMORY):       4097kB
    (hint: When tracing set SCOREP_TOTAL_MEMORY=4097kB to avoid intermediate flushes
    or reduce requirements using USR regions filters.)

    flt type max_buf[B] visits time[s] time[%] time/visit[us]  region
        ALL    377,009 14,427   65.10   100.0        4512.48  ALL
        USR    281,640 11,735    0.01     0.0           0.87  USR
        OMP     75,144  1,850    7.01    10.8        3787.55  OMP
        COM     20,184    841    0.00     0.0           4.74  COM
      SCOREP         41      1   58.08    89.2    58080341.04  SCOREP

        USR     39,240  1,635    0.00     0.0           0.28  Json::ValueType?Json::Value::type()?const
        ...
        COM      4,032    168    0.00     0.0           0.63  void?mio::abm::Simulation::evolve_world(mio::abm::TimePoint)
        COM      4,032    168    0.00     0.0           2.30  void?mio::abm::World::evolve(mio::abm::TimePoint,?mio::abm::TimeSpan)
        COM      4,032    168    0.00     0.0           6.92  void?mio::abm::World::begin_step(mio::abm::TimePoint,?mio::abm::TimeSpan)
        COM      4,032    168    0.00     0.0           5.40  void?mio::abm::World::interaction(mio::abm::TimePoint,?mio::abm::TimeSpan)
        COM      4,032    168    0.00     0.0           5.10  void?mio::abm::World::migration(mio::abm::TimePoint,?mio::abm::TimeSpan)

The output is much shorter and most functions are filtered.
If the estimated memory requirements exceed 4 GB you need to remove the responsible function from the included functions in the filter. If you added more than one function to the filter, you can check which function is responsible for the increase of required space with the output of ```scorep-score```. The functions at the top of the table are most likely to be responsible. With ```scorep-score -f filter -r scorep-folder/profile.cubex``` you can see the effect of the adjusted filter on the profile without running the experiment.
However, there are some functions that still appear in the profile that do not add overhead but are not of interest. If you checked that the overhead is acceptable you can set the filter for runtime with ```export SCOREP_FILTERING_FILE=${Path to filter}/scorep-filter-abm``` to get a clearer profile:

    Estimated aggregate size of event trace:                   94kB
    Estimated requirements for largest trace buffer (max_buf): 94kB
    Estimated memory requirements (SCOREP_TOTAL_MEMORY):       4097kB
    (hint: When tracing set SCOREP_TOTAL_MEMORY=4097kB to avoid intermediate flushes
    or reduce requirements using USR regions filters.)

    flt     type max_buf[B] visits time[s] time[%] time/visit[us]  region
        ALL     95,369  2,692   64.79   100.0       24066.57  ALL
        OMP     75,144  1,850    6.95    10.7        3758.25  OMP
        COM     20,184    841    0.00     0.0           4.98  COM
        SCOREP         41      1   57.83    89.3    57830258.74  SCOREP

         OMP     14,280    168    0.00     0.0           0.52  ?$omp?parallel?@world.cpp:68
         OMP     14,280    168    0.00     0.0           0.51  ?$omp?parallel?@world.cpp:78
         OMP     14,280    168    0.00     0.0           0.56  ?$omp?parallel?@world.cpp:152
         OMP      4,056    169    1.78     2.8       10549.08  ?$omp?for?@common_abm_loggers.h:175
         OMP      4,056    169    0.00     0.0           1.30  ?$omp?implicit?barrier?@common_abm_loggers.h:180
         OMP      4,032    168    1.54     2.4        9152.87  ?$omp?for?@world.cpp:68
         OMP      4,032    168    0.00     0.0           2.76  ?$omp?implicit?barrier?@world.cpp:73
         OMP      4,032    168    1.61     2.5        9571.18  ?$omp?for?@world.cpp:78
         OMP      4,032    168    0.00     0.0           2.75  ?$omp?implicit?barrier?@world.cpp:122
         OMP      4,032    168    2.02     3.1       12038.27  ?$omp?for?@world.cpp:152
         OMP      4,032    168    0.00     0.0           2.89  ?$omp?implicit?barrier?@world.cpp:156
         COM      4,032    168    0.00     0.0           0.70  void?mio::abm::Simulation::evolve_world(mio::abm::TimePoint)
         COM      4,032    168    0.00     0.0           2.66  void?mio::abm::World::evolve(mio::abm::TimePoint,?mio::abm::TimeSpan)
         COM      4,032    168    0.00     0.0           7.73  void?mio::abm::World::begin_step(mio::abm::TimePoint,?mio::abm::TimeSpan)
         COM      4,032    168    0.00     0.0           5.49  void?mio::abm::World::interaction(mio::abm::TimePoint,?mio::abm::TimeSpan)
         COM      4,032    168    0.00     0.0           5.17  void?mio::abm::World::migration(mio::abm::TimePoint,?mio::abm::TimeSpan)
      SCOREP         41      1   57.83    89.3    57830258.74  abm_simulation
         COM         24      1    0.00     0.0         532.58  void?mio::abm::Simulation::advance(mio::abm::TimePoint,?History&?...)??with?History?=?{mio::History<mio::abm::TimeSeriesWriter,?mio::abm::LogInfectionState>}?

### Tracing

For deeper analysis, use Scalasca for tracing:

```
export OMP_NUM_THREADS=8 #or whatever
scalasca -analyze -q -t ./bin/abm_simulation 1 <result dir> #run the program and create trace files
scalasca -examine scorep_abm_simulation_Ox8_trace #open cube to visualize the traces, this time with better information about load imbalances etc.
```

In comparison to the profiling above, it gives a timeline of events which is particularly useful for parallel programs.
Unfortunately, there is no open source tool to visualize the traces, but apparently there is a plugin for cube available (see [here](https://pramodkumbhar.com/2019/06/visualisation-of-otf2-traces-with-cubes-blade-plugin/)) and the profile becomes more accurate.
