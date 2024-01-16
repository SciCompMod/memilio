# Performance measuring with Score-P

Once you have Score-P you can build the abm with: 
```
export SCOREP_WRAPPER=off #disable scorep during cmake configure
cmake .. -DMEMILIO_ENABLE_OPENMP=ON -DCMAKE_CXX_COMPILER=${Path to scorep}/scorep-mpic++ -DCMAKE_C_COMPILER=${Path to scorep}/scorep-mpicc -DCMAKE_BUILD_TYPE=RelWithDebInfo
export SCOREP_WRAPPER=on #reenable scorep
export SCOREP_WRAPPER_INSTRUMENTER_FLAGS="--pomp --thread=omp --mpp=none --instrument-filter=${Path to filter}/filter-main"
cmake --build . -j --target abm_simulation
```

The [filter file](./filter-main) is designed to include only the interesting information about the abm. It filters (i.e. does not measure) the performance of the initialization part and some small and often called functions.

To avoid the build errors about "include style" and "gcc extension" comment out the lines that enable compiler warnings. The errors are enabled by the ```-pedantic``` flag.
It is useful to use ```--target``` to set the target to build, because Score-P has difficulties with some of the other targets.

## Sampling

To run the basic analysis with sampling:
```
export OMP_NUM_THREADS=8 #or whatever
./bin/abm_simulation 1 <some result dir>
```

There will be a directory named scorep_... with a .cubex file. That can be opened in the Cube visualization tool. This can give a first indication for bottlenecks, then you don't need Scalasca.
When changing the ABM one should check using ```scorep-score -r scorep-folder/profile.cubex``` how much memory is used by scorep and that it does not exceed 4 GB which is the maximum value for SCOREP_TOTAL_MEMORY. Also when tracing SCOREP_TOTAL_MEMORY should be adjusted to the value proposed by scorep-score. 
If the estimated memory requirements exceed 4 GB and you added small functions like getters you should add them to the excluded functions in the filter. With ```scorep-score -f filter -r scorep-folder/profile.cubex``` you can see the effect of the adjusted filter on the profile.

## Tracing

For deeper analysis, use Scalasca for tracing:

```
export OMP_NUM_THREADS=8 #or whatever
scalasca -analyze -q -t ./bin/abm_simulation 1 <result dir> #run the program and create trace files
scalasca -examine scorep_abm_simulation_Ox8_trace #open cube to visualize the traces, this time with better information about load imbalances etc.
```
