#!/bin/bash

cd ../build
## Write all data of the fictional non age resolved scenario.
dir="../data/simulation_lct_noage"
if [ ! -d "$dir" ]; then
    mkdir "$dir"
fi

# Check if all relevant folders are present.
subdir_dropR0="$dir/dropR0/"
if [ ! -d "$subdir_dropR0" ]; then
    mkdir "$subdir_dropR0"
fi
subdir_riseR0short="$dir/riseR0short/"
if [ ! -d "$subdir_riseR0short" ]; then
    mkdir "$subdir_riseR0short"
fi
subdir_riseR0shortswapped="$dir/riseR0shortswappedTETC/"
if [ ! -d "$subdir_riseR0shortswapped" ]; then
    mkdir "$subdir_riseR0shortswapped"
fi
subdir_riseR0long="$dir/riseR0long/"
if [ ! -d "$subdir_riseR0long" ]; then
    mkdir "$subdir_riseR0long"
fi

# Compile with the different numbers of subcompartments and run with different setups.
for num_subcomp in 1 3 10 50
do
    cmake -DNUM_SUBCOMPARTMENTS=$num_subcomp -DCMAKE_BUILD_TYPE="Release" .
    cmake --build . --target lct_fictional_noage

    # First case: drop R0.
    R0=0.5
    simulation_days=12
    ./bin/lct_fictional_noage $R0 $simulation_days $subdir_dropR0

    # Second case: rise R0 short 2.
    R0=2.
    simulation_days=12
    if [ "$num_subcomp" -eq 10 ] || [ "$num_subcomp" -eq 50 ]; then
        ./bin/lct_fictional_noage $R0 $simulation_days $subdir_riseR0short 1
    fi
    ./bin/lct_fictional_noage $R0 $simulation_days $subdir_riseR0short

    # Third case: rise R0 short 2 with swappedvalues for TE and TC.
    R0=2.
    simulation_days=12
    if [ "$num_subcomp" -eq 50 ]; then
        ./bin/lct_fictional_noage $R0 $simulation_days $subdir_riseR0shortswapped 1 1
    fi
    ./bin/lct_fictional_noage $R0 $simulation_days $subdir_riseR0shortswapped 0 1

    # Fourth case: Print final sizes.
    simulation_days=500
    ./bin/lct_fictional_noage 2 $simulation_days "" 0 0 1
    ./bin/lct_fictional_noage 4 $simulation_days "" 0 0 1
    ./bin/lct_fictional_noage 10 $simulation_days "" 0 0 1
done

# Fourth case: rise R0 to different R0 values and simulate for a long time period.
simulationdays=(140 100 90 80 60 60 60 60 60)
R0s=(2 3 4 5 6 7 8 9 10)
for num_subcomp in 1 2 3 4 5 6 7 8 9 10 50
do
    cmake -DNUM_SUBCOMPARTMENTS=$num_subcomp -DCMAKE_BUILD_TYPE="Release" .
    cmake --build . --target lct_fictional_noage
    for index in {0..8}
    do
        ./bin/lct_fictional_noage ${R0s[index]} ${simulationdays[index]} $subdir_riseR0long
    done
done

