#!/bin/bash

## This script can be used to get all simulation data for the numerical experiments regarding the impact of the 
## distribution assumption and of the age resolution.
## The files lct_impact_distribution_assumption.cpp and lct_impact_age_resolution.cpp are used.

# Define and construct relevant folders.
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

subdir_dropReff="$dir/dropReff/"
if [ ! -d "$subdir_dropReff" ]; then
    mkdir "$subdir_dropReff"
fi
subdir_riseReffTo2short="$dir/riseReffTo2short/"
if [ ! -d "$subdir_riseReffTo2short" ]; then
    mkdir "$subdir_riseReffTo2short"
fi
subdir_riseReffTo2shortTEhalved="$dir/riseReffTo2shortTEhalved/"
if [ ! -d "$subdir_riseReffTo2shortTEhalved" ]; then
    mkdir "$subdir_riseReffTo2shortTEhalved"
fi

# Compile with the different numbers of subcompartments and run with different setups.
for num_subcomp in 1 3 10 50
do
    cmake -DNUM_SUBCOMPARTMENTS=$num_subcomp -DCMAKE_BUILD_TYPE="Release" .
    cmake --build . --target lct_impact_distribution_assumption

    # First case: Decrease the effective reproduction number at simulation day 2 to 0.5 and simulate for 12 days.
    Reff2=0.5
    simulation_days=12
    ./bin/lct_impact_distribution_assumption $Reff2 $simulation_days $subdir_dropReff

    # Second case: Increase the effective reproduction number at simulation day 2 to 2 and simulate for 12 days.
    Reff2=2.
    simulation_days=12
    # Additionally save result with division in subcompartments.
    if [ "$num_subcomp" -eq 10 ] || [ "$num_subcomp" -eq 50 ]; then
        ./bin/lct_impact_distribution_assumption $Reff2 $simulation_days $subdir_riseReffTo2short 1
    fi
    ./bin/lct_impact_distribution_assumption $Reff2 $simulation_days $subdir_riseReffTo2short

    # Third case: Second case but with TimeExposed scaled by 0.5.
    scale_TimeExposed=0.5
    # Additionally save result with division in subcompartments.
    if [ "$num_subcomp" -eq 50 ]; then
        ./bin/lct_impact_distribution_assumption $Reff2 $simulation_days $subdir_riseReffTo2shortTEhalved 1 $scale_TimeExposed
    fi
    ./bin/lct_impact_distribution_assumption $Reff2 $simulation_days $subdir_riseReffTo2shortTEhalved 0 $scale_TimeExposed

    # Fourth case: Print final sizes without saving results.
    simulation_days=500
    ./bin/lct_impact_distribution_assumption 2 $simulation_days "" 0 1.0 1
    ./bin/lct_impact_distribution_assumption 4 $simulation_days "" 0 1.0 1
    ./bin/lct_impact_distribution_assumption 10 $simulation_days "" 0 1.0 1
done


# Fifth case: Increase the effective reproduction number at simulation day 2 to different values and 
# simulate for 200 days to compare epidemic peaks.
# Also perform simulations with TimeExposed scaled by 0.5 or 2.
# Define and construct relevant folders.
subdir_riseRefflong="$dir/riseRefflong/"
if [ ! -d "$subdir_riseRefflong" ]; then
    mkdir "$subdir_riseRefflong"
fi
subdir_riseRefflongTEhalved="$dir/riseRefflongTEhalved/"
if [ ! -d "$subdir_riseRefflongTEhalved" ]; then
    mkdir "$subdir_riseRefflongTEhalved"
fi
subdir_riseRefflongTEdoubled="$dir/riseRefflongTEdoubled/"
if [ ! -d "$subdir_riseRefflongTEdoubled" ]; then
    mkdir "$subdir_riseRefflongTEdoubled"
fi
simulationdays=200
Reff2s=(2 3 4 5 6 7 8 9 10)
for num_subcomp in 1 2 3 4 5 6 7 8 9 10 50
do
    cmake -DNUM_SUBCOMPARTMENTS=$num_subcomp -DCMAKE_BUILD_TYPE="Release" .
    cmake --build . --target lct_impact_distribution_assumption
    for index in {0..8}
    do
        ./bin/lct_impact_distribution_assumption ${Reff2s[index]} ${simulationdays} $subdir_riseRefflong
        ./bin/lct_impact_distribution_assumption ${Reff2s[index]} ${simulationdays} $subdir_riseRefflongTEhalved 0 0.5
        ./bin/lct_impact_distribution_assumption ${Reff2s[index]} ${simulationdays} $subdir_riseRefflongTEdoubled 0 2.0
    done
done

# Sixth case: Simulation for the impact of age resolution.
subdir_age_resolution="$dir/age_resolution/"
if [ ! -d "$subdir_age_resolution" ]; then
    mkdir "$subdir_age_resolution"
fi
cmake --build . --target lct_impact_age_resolution 
contact_data_dir="../../data/contacts/"
./bin/lct_impact_age_resolution $contact_data_dir $subdir_age_resolution 
