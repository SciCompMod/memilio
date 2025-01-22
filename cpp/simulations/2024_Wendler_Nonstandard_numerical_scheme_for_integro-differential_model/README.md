# A nonstandard numerical scheme for a novel SECIR integro-differential equation-based model allowing nonexponentially distributed stay times #

In this directory you find all files related to the paper
- _Anna Wendler, Lena Plötzke, Hannah Tritzschak, Martin J. Kühn (2024). A nonstandard numerical scheme for a novel SECIR integro-differential equation-based model allowing nonexponentially distributed stay times._ 
https://doi.org/10.48550/arXiv.2412.09140.

Below is an overview of the files and the paper sections they belong to. 

- Section 5.1: With the file [ide_convergence_rate](ide_convergence_rate.cpp) simulations with an integro-differential-based model (IDE model) and an ordinary differential-based model (ODE model) with different time step sizes can be performed to analyze the order of convergence of the introduced nonstandard numerical scheme. The visualization of the results can be done with [plot_convergence](plot_convergence.py).

- Section 5.2: In [ide_changepoints](ide_chanepoints.cpp) a simulation with both our IDE and an ODE model is performed where a change in the contact rate is applied after two days to investigate the model behavior at change points. The results can be visualized with [plot_changepoints](plot_changepoints.py).

- Section 5.3: The COVID-19 inspired scenario can be run with [ide_covid_inspired_scenario](ide_covid_inspired_scenario.cpp).

    For the transition probabilities we use data from [Assessment of effective mitigation and prediction of the spread of SARS-CoV-2 in Germany using demographic information and spatial resolution](https://doi.org/10.1016/j.mbs.2021.108648) and [Covasim: An agent-based model of COVID-19 dynamics and interventions](https://doi.org/10.1371/journal.pcbi.1009149). There, the data is age-resolved. We consider a model without age resolution which is why compute according parameters by averaging over the age groups. This is done in [compute_parameters](compute_parameters.cpp).

    For the stay time distributions in the respective compartments we used lognormal distributions as in [Covasim: An agent-based model of COVID-19 dynamics and interventions](https://doi.org/10.1371/journal.pcbi.1009149). There, the mean and the standard deviation were given. In order to use it with our implementation, we computed the corresponding shape and scale parameters with [get_lognormal_parameters](get_lognormal_parameters.py).

    For initialization of the models we use data reported from RKI that was averaged using a moving_average of 7.
    - Data on reported cases can be downloaded with [getCaseData](../../../pycode/memilio-epidata/memilio/epidata/getCaseData.py). To obtain the data, use the options moving_average = 7 and impute_dates = True when downloading the data. For the other function parameters, one can use default parameters from [defaultDict](../../../pycode/memilio-epidata/memilio/epidata/defaultDict.py).
    
    When simulating, we set the contact rate such that the daily new transmissions of the IDE simulation are consistent with the reported data at the beginning of the simulation. The contact scaling and the visualization is done in [run_and_plot_covid_inspired_scenario](run_and_plot_covid_inspired_scenario.py). For the scaling and the subsequent visualization we used case data from RKI (not averaged) and data on ICU patients as reported from DIVI.
    - To obtain the case data, use [getCaseData](../../../pycode/memilio-epidata/memilio/epidata/getCaseData.py) with the option impute_dates = True when downloading.
    - To obtain the data on ICU patients, use [getDIVIData](../../../pycode/memilio-epidata/memilio/epidata/getDIVIData.py) with the option impute_dates = True when downloading.

- Section 6: The share of the considered population in considered age groups that we expect to be in the ICU in June and October in comparison to the proportion of each age group in Germany is visualized using [investigate_age_distribution](investigate_age_distribution.py). This is done using age-resolved data on reported cases which can be obtained with [getCaseData](../../../pycode/memilio-epidata/memilio/epidata/getCaseData.py) with the option impute_dates = True when downloading.


# Requirements
We require the libraries `JsonCpp` and `HDF5` for running the scripts (these are optional for the whole project, see [README](../../README.md)).

The memilio.epidata package needs to be installed for the python plot scripts. Have a look at the [pycode README](../../../pycode/README.rst) and the [memilio-epidata README](../../../pycode/memilio-epidata/README.rst) for instructions how to install the package.
