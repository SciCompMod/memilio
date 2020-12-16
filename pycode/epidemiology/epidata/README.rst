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

https://opendata.arcgis.com/datasets/abad92e8eead46a4b0d252ee9438eb53_1.csv

Data from John Hopkins University (JH)

We want to get data from the Spanish Ministery of Health (MISAN) provided in the github repo:

https://github.com/datadista/datasets/tree/master/COVID%2019

Data from DIVI Intensivregister (DIVI)

Dependencies
------------

Needed python packages:

- pandas
- xlrd
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
      getdividata
      getalldata

While running the program close one figure-window to get the next one.

Run options
~~~~~~~~~~~

There are several optional run options

optional arguments working for all are:
-h, --help            show this help message and exit
  -r, --read-from-disk  Reads the data from file "json" instead of downloading
                        it.
  -ff {json,hdf5,json_timeasstring}, --file-format {json,hdf5,json_timeasstring}
                        Defines output format for data files. Default is
                        "json_timeasstring".
  -o OUT_PATH, --out-path OUT_PATH
                        Defines folder for output.

optional arguments working for some are:

  -ed END_DATE, --end-date END_DATE
                        Defines date after which data download is
                        stopped.Should have form: YYYY-mm-dd. Default is today
  -p, --plot            Plots the data.
  -sd START_DATE, --start-date START_DATE
                        Defines start date for data download. Should have
                        form: YYYY-mm-dd.Default is 2020-04-24
  -u, --update          Reads the data from file "json", downloads and adds
                        data from today.
  --split_berlin        Does not concatenate the different districts of Berlin 
			into one county and keeps it as 7 different districts
     			which are provided by the original RKI data.

Results
-------

The data is written either in json or hdf5 format

For RKI:

When speaking about infected, means always infected inclusive the already recovered persons

 ======== ======== ======================== =================
 Source   Folder   Files                    Data descritpion
 ======== ======== ======================== =================
 RKI      Germany  all_germany_rki          infected, deaths, recovered over time for whole Germany
 RKI      Germany  infected_rki             Numbers of infected over time for whole Germany
 RKI      Germany  deaths_rki               Numbers of deaths over time for whole Germany
 RKI      Germany  infected_state_rki       infected over time for different states (Bundesländer)
 RKI      Germany  all_state_rki            infected, deaths, recovered over time for different states (Bundesländer)
 RKI      Germany  infected_county_rki      infected over time for different counties (Landkreise)
 RKI      Germany  all_county_rki           infected, deaths, recovered over time for different counties (Landkreise)
 RKI      Germany  all_gender_rki           infected, deaths, recovered over time for different gender
 RKI      Germany  all_age_rki              infected, deaths, recovered over time for different age ranges
 RKI      Germany  all_state_age_rki        infected, deaths, recovered over time for different age ranges and states
 RKI      Germany  all_state_gender_rki     infected, deaths, recovered over time for different genders and states
 RKI      Germany  all_county_age_rki       infected, deaths, recovered over time for different age ranges and counties
 RKI      Germany  all_county_gender_rki    infected, deaths, recovered over time for different genders counties

 RKI-Estimation      Germany  all_germany_rki_estimated          infected, deaths, recovered, recovered_estimated, deaths_estimated over time for whole Germany
 RKI-Estimation      Germany  all_state_rki_estimated            infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different states (Bundesländer)
 RKI-Estimation      Germany  all_county_rki_estimated           infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different counties (Landkreise)
 RKI-Estimation      Germany  all_gender_rki_estimated           infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different gender
 RKI-Estimation      Germany  all_age_rki_estimated              infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different age ranges
 RKI-Estimation      Germany  all_state_age_rki_estimated        infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different age ranges and states
 RKI-Estimation      Germany  all_state_gender_rki_estimated     infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different genders and states
 RKI-Estimation      Germany  all_county_age_rki_estimated       infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different age ranges and counties
 RKI-Estimation      Germany  all_county_gender_rki_estimated    infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different genders counties


 P        Germany  FullDataB                Full data for Bundesländer
 P        Germany  FullDataL                Full data for Landkreise
 P        Germany  PopulStates              Einwohnerzahl (EWZ) for all Bundesländer
 P        Germany  PopulCounties            Einwohnerzahl (EWZ) for all Landkreise (however some are missing compared to RKI data)
 P	  Germany  county_population        Einwohnerzahl for different age groups from the 2011 census
 P	  Germany  county_current_populationEinwohnerzahl for different age groups from the 2011 census, extrapolated to the current level

 JH       .        FullData_JohnHopkins     Data as downloaded from github
 JH       .        all_provincestate        Time-cumsum of confirmed, recovered, death for states or provinces if they where given
 JH       .        all_countries            Time-cumsum of confirmed, recovered, death for every country
 JH       Germany  whole_country_Germany_jh Time-cumsum of confirmed, recovered, death for Germany
 JH       Spain    whole_country_Spain_jh   Time-cumsum of confirmed, recovered, death for Spain
 JH       France   whole_country_France_jh  Time-cumsum of confirmed, recovered, death for France
 JH       China    whole_country_China_jh   Time-cumsum of confirmed, recovered, death for China

 MISAN    Spain    spain_all_age            ['Date', 'Age', 'Gender', 'Confirmed', 'Hospitalized', 'ICU', 'Deaths'] for different age ranges
 MISAN    Spain    spain_all_state          ['Date', 'ID_State', 'State', 'Confirmed_total', 'Confirmed_PCR', 'Confirmed_AB', 'Hospitalized', 'ICU', 'Deaths', 'Recovered']
 
 DIVI     Germany  FullData_DIVI            Full data as downloaded from archive with columns ['County', 'State', 'anzahl_meldebereiche', 'reporting_hospitals', 'occupied_ICU', 'free_ICU', 'ID_State', 'Date', 'ICU', 'ICU_ventilated', 'faelle_covid_aktuell_im_bundesland', 'ID_County']
 DIVI     Germany  county_divi              ICU, ICU_ventilated over time for different counties (Landkreise) with columns ['County', 'ID_County', 'ICU', 'ICU_ventilated', 'Date']
 DIVI     Germany  state_divi               ICU, ICU_ventilated over time for different states (Bundesländer) with columns ['Date', 'ICU', 'ICU_ventilated', 'ID_State', 'State']
 DIVI     Germany  germany_divi             ICU, ICU_ventilated over time for whole Germany with columns ['Date', 'ICU', 'ICU_ventilated']
 ======== ======== ======================== =================

Some more notes
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

For DIVI:

For everyday there is one file, from which we extract the date.
However, in the beginning the data was different to the later ones.
For the first two dates, 24.4. and 25.4., there is no data for ICU_ventilated (faelle_covid_aktuell_beatmet).
For the 24.4. even has the ICU data only for each state (faelle_covid_aktuell_im_bundesland) but not for every county.
Thus, it is not yet considered in the summarized data for counties, states and whole Germany. (There are
zero entries for these dates).
Not every hospital is reporting the number of corona patients in intensive care units (ICU). The number of
reporting hospitals differs from day to day and is given in FullData_DIVI.

Notes for developers
--------------------

We use dictionaries to change the columns name to have all the names the same and are able to easily change them
If data from with other languages are used please add the dictionary in "defaultDict.py" and use the exsting one.

Note: You should not use the possibilities of pandas the access the columne with dataframe.column but instead use
datafram[column] and use th dictionaries for the column-name.

When a new script to download data is added please add the functionality to the dictionary cli_dict in the cli function in getDataIntoPandasDataFrame.py
by adding a name for it a key and adding a list with in the form [comment to print, list of used parser arguments]

If a new parser-argument has to be added, you need to add two if-loops for it to the ci-function in getDataIntoPandasDataFrame.py:
first make the parser.add_argument(...) and second to append the arg-list.

