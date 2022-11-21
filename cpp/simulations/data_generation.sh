#! /bin/bash
# Run all scripts in memilio-epidata which are necessary to run the simulation files.

# activate virutal_env 
# The path to the virtual environment has to be adjusted individually
path_virtual_env="$HOME/Documents/virtualenv/bin/activate"
if ! [[ -f $path_virtual_env ]]
then echo "Set path_virtual_env to an existing python venv." && exit
else source $path_virtual_env
fi
# path to memilio dir
cd ../..
data_dir=$PWD/data/pydata

# download data
cd "pycode/memilio-epidata"
python setup.py install
python memilio/epidata/getPopulationData.py -o $data_dir
python memilio/epidata/getDIVIData.py -o $data_dir -m 7
python memilio/epidata/getSimulationData.py -o $data_dir -m 7
python memilio/epidata/getCaseData.py -o $data_dir -m 7

echo "Generation was succesful." 