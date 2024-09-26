MEmilio Epidata Package
=======================

Content
-------

- Introduction
- Installation
- Dependencies
- Running the scripts
- Testing and Coverage
- Pylint
- Additional tools
- Notes for developers (!)
- Troubleshooting

Introduction
------------

This package provides modules and scripts to download epidemiological data from various official and inofficial sources.
A more detailed description of the sources can be found in the `epidata subfolder <memilio/epidata/README.rst>`_.

Installation
------------

Use the provided ``setup.py`` script to install the package and its dependencies.

To install the package, use (from the directory that contains ``setup.py``)

.. code:: sh

    pip install .

This copies everything required to your site-packages.

For developement of code use the command 

.. code:: sh

    pip install -e .[dev]

This command allows you to work on the code without having to reinstall the package after a change. It also installs all additional dependencies required for development and maintenance.

Dependencies
------------

Required python packages:

- pandas>=2.0.0
- matplotlib
- tables
- numpy>=1.22,<1.25
- pyarrow
- openpyxl
- xlrd
- requests
- pyxlsb
- wget

Running the scripts
-------------------

After installation the scripts can be run via the following entry points.
  - getcasedata (get case data from rki, see Results: RKI-C)
  - getpopuldata (get population data, see Results: P)
  - getjhdata (get case data from john hopkins university, see Results: JH)
  - getdividata (get ICU data from DIVI, see Results: DIVI)
  - getsimdata (get simulation data including case and vaccination data from rki, population data and ICU data, see Results: RKI-C, RKI-V, P, DIVI)
  - cleandata (deletes written files) 
  - getcommutermobility (get data about commuter mobility, see Results: BAA)
  - gettestingdata (get data about number of tests, see Results: RKI-T)
  - gethospitalizationdata (get hospitalization data from RKI, see Results: RKI-H)

For a detailed description of the run options and the resulting data files written
see the `epidata subfolder <memilio/epidata/README.rst>`_.

Testing and Coverage
--------------------

The following packages are used by the tests:

- pyfakefs (creates fake directory to test that expected folders are created and data is written)
- coverage

See Installation on how to install all these dependencies automatically.

To run the tests make 

.. code:: sh

    python -m unittest

To get the coverage report do

.. code:: sh

    python -m coverage run -m unittest
    python -m coverage report
    python -m coverage xml -o coverage_python.xml
    python -m coverage html -d coverage_python

Coverage report for actual master:

:Coverage Report: https://scicompmod.github.io/memilio/coverage/python/

Inspection via pylint
---------------------
The following packages have to be installed to run pylint:

* pylint
* pylint-json2html

See Installation on how to install all these dependencies automatically.

Run pylint with the commands

.. code:: sh
    
    memiliopylint
    pylint-json2html -f jsonextended -o build_pylint/pylint.html < build_pylint/pylint_extended.json

Pylint report for actual master:

:Pylint Report: https://dlr-sc.github.io/memilio/pylint/

Additional Tools
----------------

Some additional tools for processing or analysing data can be found in the `tools directory <tools/README.md>`_.

Notes for developers
--------------------

If a new functionality shall be added please stick to the following instructions:

When you start creating a new script:

- have a look into getDataIntoPandasDataFrame.py there the main functionality which should be used is implemented.
    - get_file is used to read in data.
    - the Conf class sets relevant download options.
    - use write_dataframe to write the pandas dataframe to file.
    - use check_dir if you want to create a new folder to write data to
- use the dictionaries in defaultDict.py to rename the existing columns of your data
    - add new column names to one of the existing language dictionaries; english, german and spanish translation exists at the moment.
    - for non-english languages always use the EngEng dictionary as the key, thus we can easily change names with just changing one line.
    - in defaultDict.py a dictionary with id, state and county name, respectively exists. Please use it.
- After renaming columns, you should not use pandas dataframe.column but instead use
  dataframe[column] where column is given by the dictionaries in defaultDict.py.
  Example: ID_County = dd.GerEng['IdLandkreis'] or dd.EngEng['idCounty'].
- For extensive operations use the progress indicator to give feedback for the user
- ALWAYS use Copy-on-Write for pandas DataFrames.
- use doxygen like comments in code as
    - add description in the beginning of the file
        - ## Header
        - # @brief name descr
        - # longer description
    - add description in the beginning of every function directly after the definition
        - start and end with """
        - add a short description to first line
        - afterwards add a longer description
        - # @param name of parameter
        - # @return type description

When you add a new script

- add a executable to the setup.py in "pycode/memilio-epidata"
- add it to the cli_dict in getDataIntoPandasDataFrame.py
    - add a meaningfull key for the new script
    - as the value add a list in the form [comment to print when script is started, list of used parser arguments (optional)]
    - if more than the default parser should be added, add these parser to the  list of used parser
- add tests
- add an entry "executablename -h" to the .github/test-py/action.yml
- add an entry "executablename -o data_dl" to the .github/workflows/main.yml
- add generated data to cleanData

Adding a new parser:

- add default value to defaultDict in defaultDict.py
- add to cli_dict in getDataIntoPandasDataFrame.py which scripts use this parser
- add an if 'new parser' in what_list and add parser.add_argument()

General
- Always add unittests
- Check test coverage report, if every new feature is covered.
- Check the pylint report just comments with "refactor" are allowed.

Troubleshooting
---------------

- HDF5 errors during installation (mostly on Windows): one of the dependencies of the epidata package requires HDF5 to be installed on the system. If HDF5 is not discovered properly, this `stack overflow thread <https://stackoverflow.com/a/67765023/1151582>`_ may help resolve the issue.
