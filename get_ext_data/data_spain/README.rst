Get github data
============

Introduction
------------

We want to get data from Spain provided in the github repo:

https://github.com/datadista/datasets/tree/master/COVID%2019

The provided data is csv.

The here provided python script (without plotting just with outputting) could be used as a cronjob to get data from Spanish Ministery of Health and get every statistic we need.

Dependencies
------------

Needed python packages:

- pandas
- tables


Running the Program
-------------------

First of all create an virtual environment by calling.

.. code:: sh

   source setup_venv.sh

The needed packages will be installed


Afterwards for getting RKI data the program can be executed, by 

.. code:: sh

   python3 getSpainDataIntoPandasFrame.py

then data is downloaded and stored in two json-files raw_*.json.
Different statistics are outputted to json or hdf5 files.


Run options
~~~~~~~~~~~

There are several optional run options:

READ_DATA [default = False]: Defines if data is downloaded (False) or just the written json-file with the full data is read and used (True).
MAKE_PLOT [default = True]: Defines if plots are shown (True) or not (False). NOTE: Plots are not implemented ! 
OUT_FORM [default = json]: Defines output format of files either json or hdf5

As an example

.. code:: sh

   python3 getSpainDataIntoPandasFrame.py READ_DATA=True



Results
-------

Following data is written either in json or hdf5 format

For Spain:

When speaking about infected, means always infected inclusive the already recovered persons

There are different columns of infected:

'Confirmed_PCR' means that these infected people were tested and confirmed to be infected by a PCR test
'Confirmed_AB' means that these infected people were tested and confirmed to be infected by an ANTIBODY test
'Confirmed_total' is the sum of the previous two
'Confirmed' if the differentiation between PCR and ANTIBODY is not made/known, only the column 'Confirmed' appears

 ======================= ================= 
 Files                   Data descritpion 
 ======================= =================
 spain_all_age             ['Date', 'Age', 'Gender', 'Confirmed', 'Hospitalized', 'ICU', 'Deaths'] for different age ranges
 spain_all_state           ['Date', 'ID_State', 'State', 'Confirmed_total', 'Confirmed_PCR', 'Confirmed_AB', 'Hospitalized', 'ICU', 'Deaths', 'Recovered']
 ======================= ================= 

 IMPORTANT NOTE: ONLY USE THIS DATA WITH CARE, WE ARE WAITING FOR AN UPDATE TO CORRECT THE FOLLOWING PROBLEM:

 ############################################################################################################
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
#                                                                                                          #
############################################################################################################