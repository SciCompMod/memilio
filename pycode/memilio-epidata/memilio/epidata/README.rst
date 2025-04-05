.. _epidata_readme:

Content
-------

- Sources
- Run options
- Results
- Notes for developers (!)

Sources
-------

- Robert Koch institute (RKI):

  - Case data (RKI-C)

    Robert Koch-Institut (2021): SARS-CoV-2 Infektionen in Deutschland, Berlin: Zenodo. DOI:10.5281/zenodo.4681153.

    We download the data from github: https://github.com/robert-koch-institut/SARS-CoV-2-Infektionen_in_Deutschland

    If the data on github is not available we download the case data from rki from
    https://npgeo-corona-npgeo-de.hub.arcgis.com/datasets/e408ccf8878541a7ab6f6077a42fd811_0


  - Vaccination data (RKI-V)

    https://github.com/robert-koch-institut/COVID-19-Impfungen_in_Deutschland

    Robert Koch-Institut (2021): COVID-19-Impfungen in Deutschland, Berlin: Zenodo. DOI:10.5281/zenodo.5126652

  - Testing Data (RKI-T)

    https://github.com/robert-koch-institut/SARS-CoV-2-PCR-Testungen_in_Deutschland/raw/main/SARS-CoV-2-PCR-Testungen_in_Deutschland.csv
  
  - Hospitalization data (RKI-H)
  
    https://raw.githubusercontent.com/robert-koch-institut/COVID-19-Hospitalisierungen_in_Deutschland/master/Aktuell_Deutschland_COVID-19-Hospitalisierungen.csv

- Population data (P) like "Einwohnerzahl" for Bundesländer and Landkreise:

  https://opendata.arcgis.com/datasets/abad92e8eead46a4b0d252ee9438eb53_1.csv

  https://www.destatis.de/DE/Themen/Laender-Regionen/Regionales/Gemeindeverzeichnis/Administrativ/04-kreise.xlsx;?__blob=publicationFile

  https://www.zensus2011.de/SharedDocs/Downloads/DE/Pressemitteilung/DemografischeGrunddaten/1A_EinwohnerzahlGeschlecht.xls?__blob=publicationFile&v=5

  https://www.regionalstatistik.de/genesis/online

- Data from DIVI Intensivregister (DIVI)

  https://www.intensivregister.de/#/aktuelle-lage/downloads

- Commuter Data from "Bundesagentur fuer Arbeit" (BAA)

  https://statistik.arbeitsagentur.de/SiteGlobals/Forms/Suche/Einzelheftsuche_Formular.html?submit=Suchen&topic_f=beschaeftigung-sozbe-krpend

- Data from "Stastistischen Bundesamt" destatis (DES)

  https://www.destatis.de/DE/Themen/Laender-Regionen/Regionales/Gemeindeverzeichnis/_inhalt.html

- Data from John Hopkins University (JH)

  https://github.com/datasets/covid-19

Running the scripts
-------------------

For informations on installation, dependencies and how to rum the scripts,
see `epidata README <../../README.rst>`_ of the above folder.

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
| -n, --no-raw                                | Defines if raw data will be stored for further use.       |
+---------------------------------------------+-----------------------------------------------------------+
| --no-progress-indicators                    | Disables all progress indicators (used for downloads etc.)|
+---------------------------------------------+-----------------------------------------------------------+
| --interactive                               | Interactive download (Handle warnings, passwords etc.).   |
+---------------------------------------------+-----------------------------------------------------------+
| -v, --verbose                              | Increases verbosity level.                                |
+---------------------------------------------+-----------------------------------------------------------+
| --skip-checks                               | Skips sanity checks etc.                                  |
+---------------------------------------------+-----------------------------------------------------------+

optional arguments working for some are:

+---------------------------------------------+-----------------------------------------------------------+
| -p, --make-plot                             | Plots the data.                                           |
+---------------------------------------------+-----------------------------------------------------------+
| -ed, --end-date                             | Changes date for which data collection is stopped         |
+---------------------------------------------+-----------------------------------------------------------+
| -sd, --start-date                           | Changes date for which data collection is started         |
+---------------------------------------------+-----------------------------------------------------------+
| -i, --impute-dates                          | Returns dataframes with all dates instead of only dates   |
|                                             | where new cases have been reported.                       |
|                                             |                                                           |
|                                             | Note that this option will have a negative impact         |
|                                             | on performance as well as on the storage space needed.    |
|                                             |                                                           |
+---------------------------------------------+-----------------------------------------------------------+
| -m N, --moving-average N                    | The central N days moving average is computed for the     |
|                                             | data.                                                     |
|                                             |                                                           |
|                                             | Note that the --impute_dates option will be implicitly    |
|                                             | turned on, as computing the moving average requires all   |
|                                             | dates to be available.                                    |
+---------------------------------------------+-----------------------------------------------------------+
| -sb, --split-berlin                         | Berlin data is split into different counties,             |
|                                             | instead of having only one county for Berlin. [cases]     |
+---------------------------------------------+-----------------------------------------------------------+
| --rep-date                                  | The reporting date will be prefered over possibly given   |
|                                             | dates of disease onset. [cases]                           |
+---------------------------------------------+-----------------------------------------------------------+
| --sanitize-data                             | Different ways to distribute vaccinations to home         |
|                                             | locations of vaccinated persons[vaccination]              |
+---------------------------------------------+-----------------------------------------------------------+
| --files                                     | Files to write [case]                                     |
+---------------------------------------------+-----------------------------------------------------------+


Hint:
When using the "--make-plot" option close one figure-window to get the next one.

Results
-------

The data is written either in json or hdf5 format

The number of "infected" persons is exported as cumulative sum such that "infected" also includes already recovered or deceased persons.
Note that for Germany, vaccinations were not reported with the home county of the vaccinated persons but with the county of vaccination.

Note for DIVI:

Not every hospital is reporting the number of corona patients in intensive care units (ICU). The number of
reporting hospitals differs from day to day and is given in FullData_DIVI.

============== ==========  =================================== =================
Source         Folder      Files                               Data description
============== ==========  =================================== =================
RKI-C          Germany     cases_infected                      numbers of infected over time for whole Germany
RKI-C          Germany     cases_deaths                        numbers of deaths over time for whole Germany
RKI-C          Germany     cases_all_germany                   infected, deaths, recovered over time for whole Germany
RKI-C          Germany     cases_infected_state                infected over time for different states (Bundesländer)
RKI-C          Germany     cases_all_state                     infected, deaths, recovered over time for different states (Bundesländer)
RKI-C          Germany     cases_infected_county               infected over time for different counties (Landkreise)
RKI-C          Germany     cases_all_county                    infected, deaths, recovered over time for different counties (Landkreise)
RKI-C          Germany     cases_all_gender                    infected, deaths, recovered over time for different gender
RKI-C          Germany     cases_all_age                       infected, deaths, recovered over time for different age ranges
RKI-C          Germany     cases_all_state_age                 infected, deaths, recovered over time for different age ranges and states
RKI-C          Germany     cases_all_state_gender              infected, deaths, recovered over time for different genders and states
RKI-C          Germany     cases_all_county_age                infected, deaths, recovered over time for different age ranges and counties
RKI-C          Germany     cases_all_county_gender             infected, deaths, recovered over time for different genders counties

RKI-V          Germany     vacc_county                         administered vaccinations per county (first, second and third shot without age resolution)
RKI-V          Germany     vacc_states                         administered vaccinations per state (first, second and third shot without age resolution)
RKI-V          Germany     vacc_county_agevacc                 administered vaccinations per county (first, second and third shot for age groups as in input
                                                               data frame, i.e., 5-11, 12-17, 18-59, 60+)
RKI-V          Germany     vacc_states_agevacc                 administered vaccinations per state (first, second and third shot for age groups as in input
                                                               data frame, i.e., 5-11, 12-17, 18-59, 60+)
RKI-V          Germany     vacc_county_ageinf                  administered vaccinations per county (first, second and third shot for age groups as in cases
                                                               data frame, i.e., 0-4, 5-14, 15-34, 35-59, 60-79, 80+)
RKI-V          Germany     vacc_states_ageinf                  administered vaccinations per state (first, second and third shot for age groups as in cases
                                                               data frame, i.e., 0-4, 5-14, 15-34, 35-59, 60-79, 80+)

RKI-T          Germany     germany_testpos                     potive rates of tests over time for germany
RKI-T          Germany     germany_states_testpos              positve rates of tests over time for different states
RKI-T          Germany     germany_conties_from_states_testpos positive rates of tests over time for different counties from positive rate for states

RKI-H          Germany     hospit_state_age                    hospitalizations per day for different age groups and states
RKI-H          Germany     hospit_germany_age                  hospitalizations per day in germany for different age groups
RKI-H          Germany     hospit_state_age                    hospitalizations per day for different states
RKI-H          Germany     hospit_germany                      hospitalizations per day in germany

P              Germany     county_current_population[_dim401]  population for different age groups from the 2011 census, extrapolated to the current level [with Wartburgkreis and Eisenach separated]
P              Germany     county_population[_dim401]          population for different age groups from the 2011 census [with Wartburgkreis and Eisenach separated]
P              Germany     county_table                        raw information on the German counties and its population sizes
P              Germany     reg_key                             unchanged regional keys from excel table
P              Germany     zensus                              unchanged zensus data

JH             .           FullData_JohnHopkins                data as downloaded from github
JH             .           all_provincestate                   time-cumsum of confirmed, recovered, death for states or provinces if they where given
JH             .           all_countries                       time-cumsum of confirmed, recovered, death for every country
JH             Germany     whole_country_Germany_jh            time-cumsum of confirmed, recovered, death for Germany
JH             Spain       whole_country_Spain_jh              time-cumsum of confirmed, recovered, death for Spain
JH             France      whole_country_France_jh             time-cumsum of confirmed, recovered, death for France
JH             Italy       whole_country_Italy_jh              time-cumsum of confirmed, recovered, death for Italy
JH             SouthKorea  whole_country_SouthKorea_jh         time-cumsum of confirmed, recovered, death for SouthKorea
JH             China       whole_country_China_jh              time-cumsum of confirmed, recovered, death for China
JH             US          whole_country_US_jh                 time-cumsum of confirmed, recovered, death for US

DIVI           Germany     FullData_DIVI                       full data as downloaded from archive with columns ['County', 'State', 'anzahl_meldebereiche', 'reporting_hospitals', 'occupied_ICU', 'free_ICU', 'ID_State', 'Date', 'ICU', 'ICU_ventilated', 'faelle_covid_aktuell_im_bundesland', 'ID_County']
DIVI           Germany     county_divi                         ICU, ICU_ventilated over time for different counties (Landkreise) with columns ['County', 'ID_County', 'ICU', 'ICU_ventilated', 'Date']
DIVI           Germany     state_divi                          ICU, ICU_ventilated over time for different states (Bundesländer) with columns ['Date', 'ICU', 'ICU_ventilated', 'ID_State', 'State']
DIVI           Germany     germany_divi                        ICU, ICU_ventilated over time for whole Germany with columns ['Date', 'ICU', 'ICU_ventilated']

BAA            Germany     mobility_bfa_2020_dim401           number of commuters from one county into another indexed by county ids (with county Eisenach)
BAA            Germany     mobility_bfa_2020_dim400           number of commuters from one county into another indexed by county ids (without county Eisenach, i.e., Eisenach merged into county Wartburgkreis)
============== ==========  =================================== =================

More detailed information can be found in the
`documentation <https://scicompmod.github.io/memilio/documentation/index.html>`_  of the different functions.

Notes for developers
--------------------

If a new functionality shall be added please stick to the instructions in `epidata README <../../README.rst>`_ of the above folder.

For information about testing, coverage, pylint and tools see also the `epidata README <../../README.rst>`_ of the above folder.
