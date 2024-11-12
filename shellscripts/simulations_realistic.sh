#!/bin/bash

cd ../build
# Define and construct all relevant folders.
data_dir="../data"
result_folder="/simulation_lct_real/"
result_dir="$data_dir$result_folder"
if [ ! -d "$result_dir" ]; then
    mkdir "$result_dir"
fi

# Define some parameters
year=2020
RelativeTransmissionNoSymptoms=1.
RiskOfInfectionFromSymptomatic=0.3
# For October:
month_oct=10
day_oct=1
scale_contacts_oct=0.6531
npi_size_oct=0.3
# For July: 
month_jul=7
day_jul=1
scale_contacts_jul=0.4103
npi_size_jul=-0.3

# Compile with the different numbers of subcompartments and run with different setups.
for num_subcomp in 1 3 10 50
do
    cmake -DNUM_SUBCOMPARTMENTS=$num_subcomp -DCMAKE_BUILD_TYPE="Release" .
    cmake --build . --target lct_realistic_scenario

    # First case: 01/10/2020.
    ./bin/lct_realistic_scenario 0 $year $month_oct $day_oct $data_dir $result_folder $RelativeTransmissionNoSymptoms $RiskOfInfectionFromSymptomatic $scale_contacts_oct $npi_size_oct

    # Second case: 01/07/2020.
    ./bin/lct_realistic_scenario 0 $year $month_jul $day_jul $data_dir $result_folder $RelativeTransmissionNoSymptoms $RiskOfInfectionFromSymptomatic $scale_contacts_jul $npi_size_jul

done

./bin/lct_realistic_scenario 1 $year $month_oct $day_oct $data_dir $result_folder $RelativeTransmissionNoSymptoms $RiskOfInfectionFromSymptomatic $scale_contacts_oct $npi_size_oct
./bin/lct_realistic_scenario 1 $year $month_jul $day_jul $data_dir $result_folder $RelativeTransmissionNoSymptoms $RiskOfInfectionFromSymptomatic $scale_contacts_jul $npi_size_jul
