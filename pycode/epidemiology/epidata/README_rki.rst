Get ArcGis data
============

Introduction
------------

We want to get data from RKI. 

RKI Dashboard: https://experience.arcgis.com/experience/478220a4c454480e823b17327b2bf1d4/page/page_1/

You can find the data also on:

https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0

The provided data is either geojson or csv.

The here provided python script (without plotting just with outputting) could be used as a cronjob to get data from rki and get every statistic we need.

We also get the "Einwoherzahl" for Bundesländer and Landkreise

https://opendata.arcgis.com/datasets/5dc2fc92850241c3be3d704aa0945d9c_2.csv

https://opendata.arcgis.com/datasets/b2e6d8854d9744ca88144d30bef06a76_1.geojson

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


Afterwards for getting RKI data the program can be executed, by 

.. code:: sh

   python3 getRKIData.py

For getting population data do

.. code:: sh

   python3 getPopulationData.py

than data is downloaded and analysed.
The full data set is stored in a json-file.
Different statistic is plotted and outputted to json or hdf5 files.


While running the program close one figure-window to get the next one.


Run options
~~~~~~~~~~~

There are several optional run options:

READ_DATA [default = False]: Defines if data is downloaded (False) or just the written json-file with the full data is read and used (True).
MAKE_PLOT [default = True]: Defines if plots are shown (True) or not (False).   
OUT_FORM [default = json]: Defines output format of files either json or hdf5

As an example

.. code:: sh

   python3 getRKIData.py READ_DATA=True



Results
-------

Following data is written either in json or hdf5 format

For RKI:

When speaking about infected, means always infected inclusive the already recovered persons


 ======================= ================= 
 Files                   Data descritpion 
 ======================= =================
 infected                Numbers of infected over "Meldedatum" for whole Germany
 deaths                  Numbers of deaths over "Meldedatum" for whole Germany
 infected_state          infected over "Meldedatum" for different states (Bundesländer)
 gbNF_state_nested       same data as above but different output
 all_state               infected, deaths, recovered over "Meldedatum" for different states (Bundesländer)
 infected_county         infected over "Meldedatum" for different counties (Landkreise)
 all_county              infected, deaths, recovered over "Meldedatum" for different counties (Landkreise)
 all_gender              infected, deaths, recovered over "Meldedatum" for different gender
 all_age                 infected, deaths, recovered over "Meldedatum" for different age ranges
 all_state_age           infected, deaths, recovered over "Meldedatum" for different age ranges and states
 all_state_gender        infected, deaths, recovered over "Meldedatum" for different genders and states
 all_county_age          infected, deaths, recovered over "Meldedatum" for different age ranges and counties
 all_county_gender       infected, deaths, recovered over "Meldedatum" for different genders counties
 ======================= ================= 


For population data
======================= ================= 
 Files                   Data descritpion
======================= =================
FullDataB               Full data for Bundesländer
FullDataL               Full data for Landkreise
PopulStates             Einwohnerzahl (EWZ) for all Bundesländer 
PopulCounties           Einwohnerzahl (EWZ) for all Landkreise (however some are missing compared to RKI data)
======================= ================= 

