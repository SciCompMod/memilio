# Revisiting the Linear Chain Trick in epidemiological models: Implications of underlying assumptions for numerical solutions #

In this directory you will find all files related to the Paper 

- _Lena Plötzke , Anna Wendler , René Schmieding , Martin J. Kühn (2024). Revisiting the Linear Chain Trick in epidemiological models: Implications of underlying assumptions for numerical solutions._ 
https://doi.org/10.48550/arXiv.2412.09140.

Below is an overview of the files and the paper sections they belong to.
- Section 4.2: The file [lct_impact_distribution_assumption](lct_impact_distribution_assumption.cpp) provides the functionality to run simulations to assess the impact of the distribution assumption or different numbers of subcompartments. The dynamics at change points and epidemic peaks can be examined using this script. The population is not divided into age groups for these experiments. All simulation results are created and saved in `.h5` files when the shell script [get_data_numerical_experiments](get_data_numerical_experiments.sh) is executed. The visualizations of these simulation results in the paper were created using the python script TODO.
- Section 4. Attention: for get_data_covid19_inspired.sh you need to download divi and rki data beforehand using epidata


For most of the above `.cpp` files, the number of subcompartments used in the LCT models for all compartments can be controlled via the preprocessor macro NUM_SUBCOMPARTMENTS. Have a look at the files for further documentation or the shell scripts for the usage. 
Attention: We require the libraries `JsonCpp` and `HDF5` for running the scripts (these are optional for the whole project, see [README](../../README.md)).
