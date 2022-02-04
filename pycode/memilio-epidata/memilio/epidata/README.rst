.. _epidata_readme:

Epidemiology python package - Epidata Subpackage
================================================

Content
-------

- Introduction
- Dependencies
- Running the scripts
- Results
- Notes for developers (!)
- Some more notes

Information
-----------

:Documentation: https://hpc-against-corona.pages.gitlab.dlr.de/epidemiology/master/documentation/index.html
:Python Coverage Report: https://hpc-against-corona.pages.gitlab.dlr.de/epidemiology/master/coverage/python/index.html
:Pylint Report: https://hpc-against-corona.pages.gitlab.dlr.de/epidemiology/master/pylint/pylint.html


Introduction
------------

Getting data from different sources and convert them to usable data using the python pandas package

Our sources are:

- Robert Koch institute (RKI) For German data:

  RKI Dashboard: https://experience.arcgis.com/experience/478220a4c454480e823b17327b2bf1d4/page/page_1/

  You can find the data also on:

  https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/dd4580c810204019a7b8eb3e0b329dd6_0

  The provided data is either geojson or csv.

- Population data (P) like "Einwoherzahl" for Bundesländer and Landkreise:

  https://opendata.arcgis.com/datasets/5dc2fc92850241c3be3d704aa0945d9c_2.csv

  https://opendata.arcgis.com/datasets/b2e6d8854d9744ca88144d30bef06a76_1.geojson

  https://www.zensus2011.de/SharedDocs/Downloads/DE/Pressemitteilung/DemografischeGrunddaten/1A_EinwohnerzahlGeschlecht.xls?__blob=publicationFile&v=5

- Data from John Hopkins University (JH)

  We want to get data from the Spanish Ministery of Health (MISAN) provided in the github repo:

  https://github.com/datadista/datasets/tree/master/COVID%2019

- Data from DIVI Intensivregister (DIVI)

- (MISAN)

Dependencies
------------

Needed python packages:

- pandas<1.2.0
- matplotlib
- tables
- numpy>=1.21
- openpyxl
- xlrd
- requests

Running the scripts
-------------------

To run the scripts use the setup.py in the folder "epidemiology/pycode/" and everything is installed and useable via several entry points.
For details see README.rst in the above folder.


Run options
~~~~~~~~~~~

There are several optional run options

optional arguments working for all are:

+---------------------------------------------+-----------------------------------------------------------+
| -h, --help                                  | show this help message and exit                           |
+---------------------------------------------+-----------------------------------------------------------+
| -r, --read-data                             | Reads the data from file "json" instead of downloading it.|
+---------------------------------------------+-----------------------------------------------------------+
| -o OUT_FOLDER,                              | Defines folder for output.                                |
| --out-folder OUT_FOLDER                     |                                                           |
+---------------------------------------------+-----------------------------------------------------------+
| -ff {json,hdf5,json_timeasstring}           | Defines output format for data files.                     |
| --file-format {json,hdf5,json_timeasstring} | Default is "json_timeasstring".                           |
+---------------------------------------------+-----------------------------------------------------------+

optional arguments working for some are:

+---------------------------------------------+-----------------------------------------------------------+
| -p, --make-plot                             | Plots the data.                                           |
+---------------------------------------------+-----------------------------------------------------------+
| -ed, --end-date                             | Changes date for which data collection is stopped [divi]  |
+---------------------------------------------+-----------------------------------------------------------+
| -sd, --start-date                           | Changes date for which data collection is started [divi]  |
+---------------------------------------------+-----------------------------------------------------------+
| -fd, --fill-dates                           | Returns dataframes with all dates instead of only dates   |
|                                             | where new cases have been reported.                       |
|                                             |  Note that this option will have a negative impact        |
|                                             |  on performance as well as on the storage space needed.   |
|                                             |  [rki]                                                    |
+---------------------------------------------+-----------------------------------------------------------+
| -ma, --moving-average                       | The 7 day moving average is computed for the data.        |
|                                             |  Note that the --fill_dates option will be implicitly     |
|                                             |  turned on, as computing the moving average requires all  |
|                                             |  dates to be available. [rki]                             |
+---------------------------------------------+-----------------------------------------------------------+
| -sb, --split-berlin                         | Berlin data is split into different counties              |
|                                             |  , instead of having only one county for Berlin. [rki]    |
+---------------------------------------------+-----------------------------------------------------------+
| -u, -- update-data                          | Just chronological missing data is added,                 |
|                                             | **after** the existing ones [divi]                        |
+---------------------------------------------+-----------------------------------------------------------+

Hint:
When using the "--make-plot" option close one figure-window to get the next one.

Results
-------

The data is written either in json or hdf5 format

When speaking about infected, means always infected inclusive the already recovered persons

 ============== ==========  ================================== =================
 Source         Folder      Files                              Data description
 ============== ==========  ================================== =================
 RKI            Germany     infected_rki                       Numbers of infected over time for whole Germany
 RKI            Germany     deaths_rki                         Numbers of deaths over time for whole Germany
 RKI            Germany     all_germany_rki                    infected, deaths, recovered over time for whole Germany
 RKI            Germany     infected_state_rki                 infected over time for different states (Bundesländer)
 RKI            Germany     all_state_rki                      infected, deaths, recovered over time for different states (Bundesländer)
 RKI            Germany     infected_county_rki                infected over time for different counties (Landkreise)
 RKI            Germany     all_county_rki                     infected, deaths, recovered over time for different counties (Landkreise)
 RKI            Germany     all_gender_rki                     infected, deaths, recovered over time for different gender
 RKI            Germany     all_age_rki                        infected, deaths, recovered over time for different age ranges
 RKI            Germany     all_state_age_rki                  infected, deaths, recovered over time for different age ranges and states
 RKI            Germany     all_state_gender_rki               infected, deaths, recovered over time for different genders and states
 RKI            Germany     all_county_age_rki                 infected, deaths, recovered over time for different age ranges and counties
 RKI            Germany     all_county_gender_rki              infected, deaths, recovered over time for different genders counties

 RKI            Germany     vaccine_data_[DATE]       administered vaccines, first shot, full vaccination, vaccination ratio, vacc ratio young, vacc ratio old

 RKI-Estimation Germany     all_germany_rki_estimated          infected, deaths, recovered, recovered_estimated, deaths_estimated over time for whole Germany
 RKI-Estimation Germany     all_state_rki_estimated            infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different states (Bundesländer)
 RKI-Estimation Germany     all_county_rki_estimated           infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different counties (Landkreise)
 RKI-Estimation Germany     all_gender_rki_estimated           infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different gender
 RKI-Estimation Germany     all_age_rki_estimated              infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different age ranges
 RKI-Estimation Germany     all_state_age_rki_estimated        infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different age ranges and states
 RKI-Estimation Germany     all_state_gender_rki_estimated     infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different genders and states
 RKI-Estimation Germany     all_county_age_rki_estimated       infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different age ranges and counties
 RKI-Estimation Germany     all_county_gender_rki_estimated    infected, deaths, recovered, recovered_estimated, deaths_estimated over time for different genders counties

 P              Germany     FullDataB                          Full data for Bundesländer
 P              Germany     FullDataL                          Full data for Landkreise
 P              Germany     PopulStates                        Einwohnerzahl (EWZ) for all Bundesländer
 P              Germany     PopulCounties                      Einwohnerzahl (EWZ) for all Landkreise (however some are missing compared to RKI data)
 P              Germany     county_population                  Einwohnerzahl for different age groups from the 2011 census
 P              Germany     county_current_population          Einwohnerzahl for different age groups from the 2011 census, extrapolated to the current level
 P              Germany     migration                          Unchanged migration data
 P              Germany     reg_key                            Unchangenged regional keys from excel table
 P              Germany     zensus                             Unchanged Zensus data

 JH             .           FullData_JohnHopkins               Data as downloaded from github
 JH             .           all_provincestate                  Time-cumsum of confirmed, recovered, death for states or provinces if they where given
 JH             .           all_countries                      Time-cumsum of confirmed, recovered, death for every country
 JH             Germany     whole_country_Germany_jh           Time-cumsum of confirmed, recovered, death for Germany
 JH             Spain       whole_country_Spain_jh             Time-cumsum of confirmed, recovered, death for Spain
 JH             France      whole_country_France_jh            Time-cumsum of confirmed, recovered, death for France
 JH             Italy       whole_country_Italy_jh             Time-cumsum of confirmed, recovered, death for Italy
 JH             SouthKorea  whole_country_SouthKorea_jh        Time-cumsum of confirmed, recovered, death for SouthKorea
 JH             China       whole_country_China_jh             Time-cumsum of confirmed, recovered, death for China
 JH             US          whole_country_US_jh                Time-cumsum of confirmed, recovered, death for US

 DIVI           Germany     FullData_DIVI                      Full data as downloaded from archive with columns ['County', 'State', 'anzahl_meldebereiche', 'reporting_hospitals', 'occupied_ICU', 'free_ICU', 'ID_State', 'Date', 'ICU', 'ICU_ventilated', 'faelle_covid_aktuell_im_bundesland', 'ID_County']
 DIVI           Germany     county_divi                        ICU, ICU_ventilated over time for different counties (Landkreise) with columns ['County', 'ID_County', 'ICU', 'ICU_ventilated', 'Date']
 DIVI           Germany     state_divi                         ICU, ICU_ventilated over time for different states (Bundesländer) with columns ['Date', 'ICU', 'ICU_ventilated', 'ID_State', 'State']
 DIVI           Germany     germany_divi                       ICU, ICU_ventilated over time for whole Germany with columns ['Date', 'ICU', 'ICU_ventilated']
 ============== ==========  ================================== =================

Notes for developers
--------------------

If a new functionality shell be added please stick to the following instructions:

When you start creating a new script:

- have a look into getDataIntoPandasDataFrame.py there the main functionality which should be used is implemented.
   - loadCsv or loadGeoJson are used to read in data
   - use the dictionaries in defaultDict.py to rename the existing columns of you data
      - add new column names to one of the existing languages; english, german and spanish translation exists at the moment.
      - for non-english languages always use the EngEng dictionary as the key, thus we can easily change names with just changing one line.
      - in defaultDict.py a dictionary with id and state and county name, respectivly exists. Please use it.
- After renaming columns, you should not use the possibilities of pandas the access the column with dataframe.column but instead use
datafram[column] and use th dictionaries to define the column-name. Example: Altersgruppe2 = dd.GerEng['Altersgruppe2']; again in this way it is easier to change the column names.
- use check_dir of getDataIntoPandasDataFrame.py if you want to create a new folder to write data to
- use write_dataframe of getDataIntoPandasDataFrame.py to write the pandas dataframe to file.
- use doxygen like comments in code as
    - add description in the beginning of the file
        - ## Header
        - # @brief name descr
        - # longer description
    - add description in the beginning of every function directly after the definiton
        - start and end with """
        - add a short description to first line
        - afterwards add a longer description
        - # @param name of parameter
        - # @return type description

When you add a new script

- add a executable to the setup.py in "epidemiology/pycode/"
- add it to the cli_dict in getDataIntoPandasDataFrame.py
    - add a meaningfull key for the new script
    - as the value add a list in the form [comment to print when script is started, list of used parser arguments (optional)]
    - if more than the default parser should be added, add these parser to the  list of used parser
- add tests
- add an entry "executablename -h" to the .gitlab-ci.yml
- add it to getAll.py
- add generated data to cleanData

Adding a new parser:

- add default value to defaultDict in defaultDict.py
- add to cli_dict in getDataIntoPandasDataFrame.py which scripts use this parser
- add an if 'new parser' in what_list and add parser.add_argument()

General
- Always add unittests
- Check test coverage report, if every new feature is covered.
- Check the pylint report just comments with "refactor" are allowed.

More detailed information can be found in the documentation of the different functions in

Some more notes
---------------

When speaking about infected, means always infected inclusive the already recovered persons

There are different columns of infected:

'Confirmed_PCR' means that these infected people were tested and confirmed to be infected by a PCR test
'Confirmed_AB' means that these infected people were tested and confirmed to be infected by an ANTIBODY test
'Confirmed_total' is the sum of the previous two
'Confirmed' if the differentiation between PCR and ANTIBODY is not made/known, only the column 'Confirmed' appears


For DIVI:

For everyday there is one file, from which we extract the date.
However, in the beginning the data was different to the later ones.
For the first two dates, 24.4. and 25.4., there is no data for ICU_ventilated (faelle_covid_aktuell_beatmet).
For the 24.4. even has the ICU data only for each state (faelle_covid_aktuell_im_bundesland) but not for every county.
Thus, it is not yet considered in the summarized data for counties, states and whole Germany. (There are
zero entries for these dates).
Not every hospital is reporting the number of corona patients in intensive care units (ICU). The number of
reporting hospitals differs from day to day and is given in FullData_DIVI.
