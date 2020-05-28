Epidemiology python package - Epidata Subpackage
================================================

Introduction
------------

Getting data from different sources and convert them to usable data

Our sources are:

From ArcGis we get the following data:

Robert Koch institute (RKI) For German data:

RKI Dashboard: https://experience.arcgis.com/experience/478220a4c454480e823b17327b2bf1d4/page/page_1/

You can find the data also on:

https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0

The provided data is either geojson or csv.

Population data like "Einwoherzahl" for Bundesländer and Landkreise:

https://opendata.arcgis.com/datasets/5dc2fc92850241c3be3d704aa0945d9c_2.csv

https://opendata.arcgis.com/datasets/b2e6d8854d9744ca88144d30bef06a76_1.geojson

Data from John Hopkins University (JH)

We want to get data from the Spanish Ministery of Health (MISAN) provided in the github repo:

https://github.com/datadista/datasets/tree/master/COVID%2019

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

The needed packages can be installed using the setup.py in the folder "/home/devel/epidemiology/pycode/"

For usage:

.. code:: sh

    python setup.py install

For developement:

.. code:: sh

    python setup.py develop

Afterwards for getting RKI data the program can be executed, by calling one of the following comments:

      getrkidata
      getpopuldata
      getjhdata
      getspaindata
      getalldata

While running the program close one figure-window to get the next one.

Run options
~~~~~~~~~~~

There are several optional run options

optional arguments:
  -h, --help                         show this help message and exit
  -r, --read-from-disk               Reads the data from file "json" instead of downloading it.
  -p, --plot                         Plots the data.
  -h5, --hdf5                        Changes output format from json to hdf5.
  -o OUT_PATH, --out_path OUT_PATH   Defines folder for output.


Note: The plot option is at the moment just working for the rki data

Results
-------

Following data is written either in json or hdf5 format

For RKI:

When speaking about infected, means always infected inclusive the already recovered persons

 ======== ======== ======================== =================
 Source   Folder   Files                    Data descritpion
 ======== ======== ======================== =================
 RKI      Germany  infected_rki             Numbers of infected over date for whole Germany
 RKI      Germany  deaths_rki               Numbers of deaths over date for whole Germany
 RKI      Germany  infected_state_rki       infected over date for different states (Bundesländer)
 RKI      Germany  all_state_rki            infected, deaths, recovered over date for different states (Bundesländer)
 RKI      Germany  infected_county_rki      infected over date for different counties (Landkreise)
 RKI      Germany  all_county_rki           infected, deaths, recovered over date for different counties (Landkreise)
 RKI      Germany  all_gender_rki           infected, deaths, recovered over date for different gender
 RKI      Germany  all_age_rki              infected, deaths, recovered over date for different age ranges
 RKI      Germany  all_state_age_rki        infected, deaths, recovered over date for different age ranges and states
 RKI      Germany  all_state_age5_rki       infected, deaths, recovered over date for different age difference of 10 years and states
 RKI      Germany  all_state_age10_rki      infected, deaths, recovered over date for different age difference of 10 and states
 RKI      Germany  all_state_gender_rki     infected, deaths, recovered over date for different genders and states
 RKI      Germany  all_county_age_rki       infected, deaths, recovered over date for different age ranges and counties
 RKI      Germany  all_county_age5_rki      infected, deaths, recovered over date for different age ranges (5 years) and counties
 RKI      Germany  all_county_age10_rki     infected, deaths, recovered over date for different age ranges (10 years) and counties
 RKI      Germany  all_county_gender_rki    infected, deaths, recovered over date for different genders counties

 P        Germany  FullDataB                Full data for Bundesländer
 P        Germany  FullDataL                Full data for Landkreise
 P        Germany  PopulStates              Einwohnerzahl (EWZ) for all Bundesländer
 P        Germany  PopulCounties            Einwohnerzahl (EWZ) for all Landkreise (however some are missing compared to RKI data)

 JH       .        FullData_JohnHopkins     Data as downloaded from github
 JH       .        all_provincestate        Time-cumsum of confirmed, recovered, death for states or provinces if they where given
 JH       .        all_countries            Time-cumsum of confirmed, recovered, death for every country
 JH       Germany  whole_country_Germany_jh Time-cumsum of confirmed, recovered, death for Germany
 JH       Spain    whole_country_Spain_jh   Time-cumsum of confirmed, recovered, death for Spain
 JH       France   whole_country_France_jh  Time-cumsum of confirmed, recovered, death for France
 JH       China    whole_country_China_jh   Time-cumsum of confirmed, recovered, death for China

 MISAN    Spain    spain_all_age            ['Date', 'Age', 'Gender', 'Confirmed', 'Hospitalized', 'ICU', 'Deaths'] for different age ranges
 MISAN    Spain    spain_all_state          ['Date', 'ID_State', 'State', 'Confirmed_total', 'Confirmed_PCR', 'Confirmed_AB', 'Hospitalized', 'ICU', 'Deaths', 'Recovered']
 ======== ======== ======================== =================

Some more Notes
---------------

When speaking about infected, means always infected inclusive the already recovered persons

There are different columns of infected:

'Confirmed_PCR' means that these infected people were tested and confirmed to be infected by a PCR test
'Confirmed_AB' means that these infected people were tested and confirmed to be infected by an ANTIBODY test
'Confirmed_total' is the sum of the previous two
'Confirmed' if the differentiation between PCR and ANTIBODY is not made/known, only the column 'Confirmed' appears

For RKI:

When the plot option is turned on: while running the program close one figure-window to get the next one.

For Spain:

IMPORTANT NOTE: ONLY USE THIS DATA WITH CARE, WE ARE WAITING FOR AN UPDATE TO CORRECT THE FOLLOWING PROBLEM:

#                                                                                                          #
#        DO NOT USE DATA FROM THE FOLLOWING REGIONS SINCE THE COLUMNS HOSPITALIZED AND ICU                 #
#        ARE NOT CORRECTLY SUMMED TO TOTAL NUMBERS ! THE SAME APPLIES TO ALL AGE DATA AT THE MOMENT !      #
#                                                                                                          #
#               HOSPITALIZED                                   ICU                                         #
#               Castilla La Mancha (until 2020-04-11)          Castilla La Mancha (hasta 2020-04-12)       #
#               Comunidad Valenciana (hasta 2020-04-08)        Castilla y León (hasta 2020-04-17)          #
#               Madrid (hasta 2020-04-26)                      Comunidad Valenciana (hasta 2020-04-08)     #
#               Castilla y León (hasta 2020-04-06)             Galicia (hasta 2020-04-29)                  #
#               Madrid (hasta 2020-04-26)                                                                  #

Notes for developers
--------------------

We use dictionaries to change the columns name to have all the names the same and are able to easily change them
If data from with other languages are used please add the dictionary in "defaultDict.py" and use the exsting one.
Note: You should not use the possibilities of pandas the access the columne with dataframe.column but instead use
datafram[column] and use th dictionaries for the column-name.



