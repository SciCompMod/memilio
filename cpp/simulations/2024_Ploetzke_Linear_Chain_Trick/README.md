# Revisiting the Linear Chain Trick in epidemiological models: Implications of underlying assumptions for numerical solutions #

In this directory you will find all files related to the Paper 

- _Lena Plötzke , Anna Wendler , René Schmieding , Martin J. Kühn (2024). Revisiting the Linear Chain Trick in epidemiological models: Implications of underlying assumptions for numerical solutions._ 
https://doi.org/10.48550/arXiv.2412.09140.

Below is an overview of the files and the paper sections they belong to.

- Section 4.2: The file [lct_impact_distribution_assumption](lct_impact_distribution_assumption.cpp) provides the functionality to run simulations to assess the impact of the distribution assumption or different numbers of subcompartments. The dynamics at change points and epidemic peaks can be examined using this script. The population is not divided into age groups for these experiments. All simulation results are created and saved in `.h5` files when the shell script [get_data_numerical_experiments](get_data_numerical_experiments.sh) is executed. The visualizations of these simulation results in the paper were created using the python script [plot_numerical_experiments](plot_numerical_experiments.py).

- Section 4.3: With the file [lct_impact_age_resolution](lct_impact_age_resolution.cpp), one can run simulations to assess the impact of including an age resolution. The simulation results are created and saved together with the results for section 4.2 with the shellscript [get_data_numerical_experiments](get_data_numerical_experiments.sh). The visualizations are also created with [plot_numerical_experiments](plot_numerical_experiments.py).

- Section 4.4: Run time measurements are possible with the file [lct_runtime](lct_runtime.cpp). `OpenMP` is used to measure the run times. The Slurm script [get_runtimes_lct](get_runtimes_lct.sh) can be used to define a job to measure the run time for different numbers of subcompartments. This script can be easily adapted to use an adaptive solver. To use the optimization flag `-O0`, uncomment the suitable line in the [CMakeLists file](CMakeLists.txt). Visualizations of the run times in the paper were created using the python script [plot_runtimes_lct](plot_runtimes_lct.py).

- Section 4.5: A COVID-19 inspired scenario in Germany in 2020 is defined in the file [lct_covid19_inspired_scenario](lct_covid19_inspired_scenario.cpp). The simulation results are created and saved with the shell script [get_data_covid19_inspired](get_data_covid19_inspired.sh). 
    The simulation is initialized using real data, which has to be downloaded beforehand:
    1. Reported case data from the RKI can be downloaded via [getCaseData](../../../pycode/memilio-epidata/memilio/epidata/getCaseData.py).

    2. Reported number of patients in intensive care units (DIVI data) can be downloaded via [getDIVIData](../../../pycode/memilio-epidata/memilio/epidata/getDIVIData.py).

    - The option impute_dates should be set to True and moving_average=0 while downloading both datasets. For the other function parameters, one can use default parameters from [defaultDict](../../../pycode/memilio-epidata/memilio/epidata/defaultDict.py) (and matching dates).

    The visualizations of the simulation results in the paper were created using the python script [plot_covid19_inspired](plot_covid19_inspired.py).

- Figure 2 and Figure 12: These figures are not based on simulation results. Figure 2 contains a visualization of the density and the survival function of Erlang distributions with different parameter choices. Figure 12 shows the age-resolved contact pattern for Germany. Both plots are created using [plot_details](plot_details.py).


For most of the above `.cpp` files, the number of subcompartments used in the LCT models for all compartments can be controlled via the preprocessor macro NUM_SUBCOMPARTMENTS. Have a look at the files for further documentation or the shell scripts for the usage. 

# Requirements
We require the libraries `JsonCpp` and `HDF5` for running the scripts (these are optional for the whole project, see [README](../../README.md)).

The memilio.epidata package needs to be installed for the python plot scripts. 
Have a look at the [pycode README](../../../pycode/README.rst) and the [memilio-epidata README](../../../pycode/memilio-epidata/README.rst) for instructions how to install the package.
