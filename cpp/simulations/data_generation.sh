#! /bin/bash
# Run all scripts in memilio-epidata which are necessary to run the simulation files.

# activate virutal_env 
# The path to the virtual environment has to be adjusted individually
path_virtual_env="$HOME/Documents/virtualenv/corona/bin/activate"
if ! [[ -f $path_virtual_env ]]
then echo "Set path_virtual_env to an existing python venv." && exit
else source $path_virtual_env
fi

# path to memilio dir
PWD=$(pwd)
cd ../..

# download data
cd "pycode/memilio-epidata"
python setup.py install
python memilio/epidata/getPopulationData.py
python memilio/epidata/getDIVIData.py -m 7
python memilio/epidata/getSimulationData.py -m 7
python memilio/epidata/getCaseData.py -m 7
python memilio/epidata/transformMobilityData.py

echo "Generation was succesful." 