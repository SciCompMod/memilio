#! /bin/bash
# Run all scripts in memilio-epidata which are necessary to run the simulation files.

# activate virtual_env 
# The path to the virtual environment activation function has to be adjusted manually
# Otherwise the script can be called with -PATH_ENV Path argument, e.g.
# sh data_generation.sh -PATH_ENV "YOUR/PATH/TO/VIRTUAL/ENV/activate"
path_virtual_env="/localdata1/kueh_mj/virtual_envs/corona-py388/bin/activate" 

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
then echo "Set path to virtual environment activation function to an existing python venv." && exit
else source $path_virtual_env
fi

# path to MEmilio dir (assumes execution of the script from the `simulations` folder)
cd ../..
data_dir=$PWD/data/pydata

# download data
cd "pycode/memilio-epidata"
python setup.py install
python memilio/epidata/getSimulationData.py -o $data_dir -m 7
python memilio/epidata/transformMobilityData.py -o $data_dir

echo "Generation was succesful." 