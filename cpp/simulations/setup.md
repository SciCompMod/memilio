# How to setup the paper simulation from brunswick

1. You need the file "braunschweig_result_ffa8.csv" downloadable a this URL:
2. You need to save this file into the folder:"/memilio/data/mobility/braunschweig_result_ffa8.csv"
3. You need to run cleanup_data.py (install numpy and pandas beforehand) please change the folders
4. You need to download the simulation files. Please follow the installation instructions at "<https://github.com/SciCompMod/memilio/tree/main/pycode/memilio-epidata>"
5. Run "python getSimulationData.py -s 2021-01-01 -e 2021-07-01 -m 7"
6. Copy the Germany folder into data/mobility/Germany
7. mkdir data/results
8. Put rki files in results
9. Change the folder in paper_abm_testing to your data folder.
10. cmake --build /Users/saschakorf/Documents/Arbeit.nosynch/memilio/memilio/cpp/build --config Release --target paper_abm_bs_testing -j 6 --
11. ./memilio/cpp/build/bin/paper_abm_bs_testing
