# Performance measuring with Score-P

Once you have Score-P you can build the abm with: 
```
export SCOREP_WRAPPER=off #disable scorep during cmake configure
cmake .. -DMEMILIO_ENABLE_OPENMP=ON -DCMAKE_CXX_COMPILER=${Path to scorep}/scorep-mpic++ -DCMAKE_C_COMPILER=${Path to scorep}/scorep-mpicc -DCMAKE_BUILD_TYPE=RelWithDebInfo
export SCOREP_WRAPPER=on #reenable scorep
export SCOREP_WRAPPER_INSTRUMENTER_FLAGS="--pomp --thread=omp --mpp=none --instrument-filter=./filter-main"
cmake --build . -j --target abm_simulation
```

For the filter file see below.

To avoid the build errors about "include style" and "gcc extension" comment out the lines that enable compiler warnings. The errors are enabled by the ```-pedantic``` flag.
It is useful to use ```--target``` to set the target to build, because Score-P has difficulties with some of the other targets.

## Sampling

To run the basic analysis with sampling:
```
export OMP_NUM_THREADS=8 #or whatever
./bin/abm_simulation 1 <some result dir>
```

There will be a directory named scorep_... with a .cubex file. That can be opened in the Cube visualization tool. This can give a first indication for bottlenecks, then you don't need Scalasca.

## Tracing

For deeper analysis, use Scalasca for tracing:
```
export OMP_NUM_THREADS=8 #or whatever
scalasca -analyze -q -t ./bin/abm_simulation 1 <result dir> #run the program and create trace files
scalasca -examine scorep_abm_simulation_Ox8_trace #open cube to visualize the traces, this time with better information about load imbalances etc.
```

## Filter file


