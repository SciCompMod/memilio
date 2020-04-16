Get RKI data
============

Introduction
------------

We want to get data from RKI. 

You can find the data also on:

https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0

The provided python script (without plotting just with outputting) could be used as a cronjob to get data from rki and get every statistic we need

Dependencies
------------

Needed python packages:
- pandas
- matplotlib


Running the Program
-------------------

First of all create an virtual environment by calling

.. code:: sh

   source setup_venv.sh


Call 

.. code:: sh

   python3 getGeoJsonIntoPandasDataFrame_for_RKI.py

than data is downloaded and analysed.
The full data set is stored in a json.
Different statistic is plotted and outputtes to json files.


While running the program close one window to get the next one.

Call

.. code:: sh

   python3 getGeoJsonIntoPandasDataFrame_for_RKI.py READ_DATA=True

To use stored data and not download it all the time.

Use 

.. code:: sh

   python3 getGeoJsonIntoPandasDataFrame_for_RKI.py MAKE_PLOT=False (READ_DATA=True)

to not see the plots just get the files.

Results
-------

When speaking about infected, means always infected inclusive the already recovered persons



Following data is written


| Files | Data descritpion |
| ------- | ------ |
| infected.json | Numbers of infected over "Meldedatum" for whole Germany |
| deaths.json | Numbers of deaths over "Meldedatum" for whole Germany |
| infected_state.json | infected over "Meldedatum" for different states (Bundesländer) |
| gbNF_state_nested.json | same data as above but different output |
| all_state.json | infected, deaths, recovered over "Meldedatum" for different states (Bundesländer)  |
| infected_county.json | infected over "Meldedatum" for different counties (Landkreise) |
| all_county.json |  infected, deaths, recovered over "Meldedatum" for different counties (Landkreise)  |
| all_gender.json |  infected, deaths, recovered over "Meldedatum" for different genders |
| all_age.json |  infected, deaths, recovered over "Meldedatum" for different age ranges |



