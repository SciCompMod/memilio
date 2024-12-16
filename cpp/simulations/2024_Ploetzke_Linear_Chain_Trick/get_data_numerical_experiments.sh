#!/bin/bash

## This script can be used to get all simulation data for the numerical experiments 
## using the files lct_impact_distribution_assumption.cpp and lct_impact_age_resolution.cpp.

# Define and construct all relevant folders.
cd ../../
if [ ! -d "build/" ]; then
    mkdir "build/"
fi
cd build/
cmake ..

dir="../../data/simulation_lct_numerical_experiments"
if [ ! -d "$dir" ]; then
    mkdir "$dir"
fi

subdir_dropR0="$dir/dropR0/"
if [ ! -d "$subdir_dropR0" ]; then
    mkdir "$subdir_dropR0"
fi
subdir_riseR0short="$dir/riseR0short/"
if [ ! -d "$subdir_riseR0short" ]; then
    mkdir "$subdir_riseR0short"
fi
subdir_riseR0shortTEhalved="$dir/riseR0shortTEhalved/"
if [ ! -d "$subdir_riseR0shortTEhalved" ]; then
    mkdir "$subdir_riseR0shortTEhalved"
fi

# Compile with the different numbers of subcompartments and run with different setups.
for num_subcomp in 1 3 10 50
do
    cmake -DNUM_SUBCOMPARTMENTS=$num_subcomp -DCMAKE_BUILD_TYPE="Release" .
    cmake --build . --target lct_impact_distribution_assumption

    # First case: drop R0.
    R0=0.5
    simulation_days=12
    ./bin/lct_impact_distribution_assumption $R0 $simulation_days $subdir_dropR0

    # Second case: rise R0 short 2.
    R0=2.
    simulation_days=12
    if [ "$num_subcomp" -eq 10 ] || [ "$num_subcomp" -eq 50 ]; then
        ./bin/lct_impact_distribution_assumption $R0 $simulation_days $subdir_riseR0short 1
    fi
    ./bin/lct_impact_distribution_assumption $R0 $simulation_days $subdir_riseR0short

    # Third case: rise R0 short 2 with TE scaled by 0.5.
    R0=2.
    simulation_days=12
    if [ "$num_subcomp" -eq 50 ]; then
        ./bin/lct_impact_distribution_assumption $R0 $simulation_days $subdir_riseR0shortTEhalved 1 0.5
    fi
    ./bin/lct_impact_distribution_assumption $R0 $simulation_days $subdir_riseR0shortTEhalved 0 0.5

    # Fourth case: Print final sizes.
    simulation_days=500
    ./bin/lct_impact_distribution_assumption 2 $simulation_days "" 0 1.0 1
    ./bin/lct_impact_distribution_assumption 4 $simulation_days "" 0 1.0 1
    ./bin/lct_impact_distribution_assumption 10 $simulation_days "" 0 1.0 1
done


# Fourth case: rise R0 to different R0 values and simulate for a long time period.
# Check if all relevant folders are present.
subdir_riseR0long="$dir/riseR0long/"
if [ ! -d "$subdir_riseR0long" ]; then
    mkdir "$subdir_riseR0long"
fi
subdir_riseR0longTEhalved="$dir/riseR0longTEhalved/"
if [ ! -d "$subdir_riseR0longTEhalved" ]; then
    mkdir "$subdir_riseR0longTEhalved"
fi
subdir_riseR0longTEdoubled="$dir/riseR0longTEdoubled/"
if [ ! -d "$subdir_riseR0longTEdoubled" ]; then
    mkdir "$subdir_riseR0longTEdoubled"
fi

simulationdays=200
R0s=(2 3 4 5 6 7 8 9 10)
for num_subcomp in 1 2 3 4 5 6 7 8 9 10 50
do
    cmake -DNUM_SUBCOMPARTMENTS=$num_subcomp -DCMAKE_BUILD_TYPE="Release" .
    cmake --build . --target lct_impact_distribution_assumption
    for index in {0..8}
    do
        ./bin/lct_impact_distribution_assumption ${R0s[index]} ${simulationdays} $subdir_riseR0long
        ./bin/lct_impact_distribution_assumption ${R0s[index]} ${simulationdays} $subdir_riseR0longTEhalved 0 0.5
        ./bin/lct_impact_distribution_assumption ${R0s[index]} ${simulationdays} $subdir_riseR0longTEdoubled 0 2.0
    done
done

subdir_age_resolution="$dir/age_resolution/"
if [ ! -d "$subdir_age_resolution" ]; then
    mkdir "$subdir_age_resolution"
fi
# Fifth case: Simualtion to check the impact of the age resolution.
cmake --build . --target lct_impact_age_resolution 
dir="../../data/contacts"
./bin/lct_impact_age_resolution $subdir_age_resolution
