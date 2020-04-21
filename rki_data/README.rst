Get RKI data
============

Introduction
------------

We want to get data from RKI. 

RKI Dashboard: https://experience.arcgis.com/experience/478220a4c454480e823b17327b2bf1d4/page/page_1/

You can find the data also on:

https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0

The provided data is either geojson or csv.

The here provided python script (without plotting just with outputting) could be used as a cronjob to get data from rki and get every statistic we need.

Dependencies
------------

Needed python packages:

- pandas
- matplotlib
- tables


Running the Program
-------------------

First of all create an virtual environment by calling.

.. code:: sh

   source setup_venv.sh

The needed packages will be installed


Afterwards the program can be executed, by do 

.. code:: sh

   python3 getGeoJsonIntoPandasDataFrame_for_RKI.py

than data is downloaded and analysed.
The full data set is stored in a json-file.
Different statistic is plotted and outputted to json files.


While running the program close one figure-window to get the next one.


Run options
~~~~~~~~~~~

There are several optional run options:

READ_DATA [default = False]: Defines if data is downloaded (False) or just the written json-file with the full data is read and used (True).
MAKE_PLOT [default = True]: Defines if plots are shown (True) or not (False).   

As an example

.. code:: sh

   python3 getGeoJsonIntoPandasDataFrame_for_RKI.py READ_DATA=True



Results
-------

When speaking about infected, means always infected inclusive the already recovered persons


Following data is written

 ======================= ================= 
 Files                   Data descritpion 
 ======================= =================
 infected.json           Numbers of infected over "Meldedatum" for whole Germany
 deaths.json             Numbers of deaths over "Meldedatum" for whole Germany
 infected_state.json     infected over "Meldedatum" for different states (Bundesländer)
 gbNF_state_nested.json  same data as above but different output
 all_state.json          infected, deaths, recovered over "Meldedatum" for different states (Bundesländer)
 all_state.h5            infected, deaths, recovered over "Meldedatum" for different states (Bundesländer), same as above but different output format 
 infected_county.json    infected over "Meldedatum" for different counties (Landkreise)
 all_county.json         infected, deaths, recovered over "Meldedatum" for different counties (Landkreise)
 all_gender.json         infected, deaths, recovered over "Meldedatum" for different gender
 all_age.json            infected, deaths, recovered over "Meldedatum" for different age ranges
 all_state_age.json      infected, deaths, recovered over "Meldedatum" for different age ranges and states
 all_state_gender.json   infected, deaths, recovered over "Meldedatum" for different genders and states
 all_county_age.json     infected, deaths, recovered over "Meldedatum" for different age ranges and counties
 all_county_gender.json  infected, deaths, recovered over "Meldedatum" for different genders counties
 ======================= ================= 



