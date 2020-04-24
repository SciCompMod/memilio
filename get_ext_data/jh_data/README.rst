Get data from John Hopkins
==========================

Introduction
------------

We want to get data from John Hopkins


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

   python3 getJHDataIntoPandasDataFrame_for_RKI.py

than data is downloaded and analysed.
The full data set is stored in a json-file.
Different statistic is plotted and outputted to json files.


While running the program close one figure-window to get the next one.


Run options
~~~~~~~~~~~

There are several optional run options:

READ_DATA [default = False]: Defines if data is downloaded (False) or just the written json-file with the full data is read and used (True).
MAKE_PLOT [default = True]: Defines if plots are shown (True) or not (False).   



Results
-------

When speaking about confirmed, means infected inclusive the already recovered persons


Following data is written

 ========================= ================= 
 Files                     Data descritpion 
 ========================= =================
 FullData_JohnHopkins.json Data as downloaded from github
 all_countries.json        Time-cumsum of confirmed, recovered, death for every country
 all_provincestate.json    Time-cumsum of confirmed, recovered, death for states or provinces if they where given
 ========================= ================= 




