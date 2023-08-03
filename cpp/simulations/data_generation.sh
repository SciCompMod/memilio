#! /bin/bash
# Run all scripts in memilio-epidata which are necessary to run the simulation files.

# activate virtual_env 
# The path to the virtual environment has to be adjusted individually
# Otherwise the script can be called with -PATH_ENV Path argument, e.g.
# sh data_generation.sh -PATH_ENV "YOUR/PATH/TO/VIRTUAL/ENV/activate"
path_virtual_env= "YOUR/PATH/TO/VIRTUAL/ENV/activate" 

# Use/check for input argument
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        -PATH_ENV)
        path_virtual_env="$2"
        shift
        shift
        ;;
        *)
        echo "Unknown parameter: $1"
        exit 1
        ;;
    esac
done

if ! [[ -f $path_virtual_env ]]
then echo "Set path_virtual_env to an existing python venv." && exit
else source $path_virtual_env
fi

# path to memilio dir
cd ../..
data_dir=$PWD/data/pydata
mobility_dir=$PWD/data/mobility/

# download data
cd "pycode/memilio-epidata"
python setup.py install
python memilio/epidata/getSimulationData.py -o $data_dir -m 7
python memilio/epidata/transformMobilityData.py -o $mobility_dir

echo "Generation was succesful." 