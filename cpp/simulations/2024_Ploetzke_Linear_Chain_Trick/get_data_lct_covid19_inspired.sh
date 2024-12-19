#!/bin/bash

## This script can be used to get all simulation data using the lct_covid19_inspired_scenario.cpp file.
## It is necessary to download RKI and DIVI data beforehand. For more information see the README.

# Define and construct all relevant folders.
cd ../../
if [ ! -d "build/" ]; then
    mkdir "build/"
fi
cd build/
data_dir="../../data"
result_folder="/simulation_lct_covid19/"
result_dir="$data_dir$result_folder"
if [ ! -d "$result_dir" ]; then
    mkdir "$result_dir"
fi
cmake ..

# Define some parameters
year=2020
RelativeTransmissionNoSymptoms=1.
RiskOfInfectionFromSymptomatic=0.3
month_oct=10
day_oct=1
scale_contacts_oct=0.6537
npi_size_oct=0.3
scale_confirmed_cases_oct=1.2

# Compile with the different numbers of subcompartments.
for num_subcomp in 1 3 10 50
do
    cmake -DNUM_SUBCOMPARTMENTS=$num_subcomp -DCMAKE_BUILD_TYPE="Release" .
    cmake --build . --target lct_covid19_inspired_scenario

    # Simulation for 01/10/2020.
    ./bin/lct_covid19_inspired_scenario 0 $year $month_oct $day_oct $data_dir $result_folder $RelativeTransmissionNoSymptoms $RiskOfInfectionFromSymptomatic $scale_contacts_oct $npi_size_oct $scale_confirmed_cases_oct

done

# Set the number of subcompartments so that each corresponds to the approximate stay time in the compartment.
./bin/lct_covid19_inspired_scenario 1 $year $month_oct $day_oct $data_dir $result_folder $RelativeTransmissionNoSymptoms $RiskOfInfectionFromSymptomatic $scale_contacts_oct $npi_size_oct $scale_confirmed_cases_oct
