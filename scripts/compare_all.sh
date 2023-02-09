#! /bin/bash
# compares the sequential and parallel output of paper_202011

tolerance=1e-10
# start date must be same as simulation start date
startdate="2020-12-12"
# paths must be absolute, or relative to the git project's root directory
data="data"
output="rounding"
saves="$output/save"
results="$output/result"

virtual_env="$HOME/Documents/virtualenv/bin/activate"

binary="build/simulations/paper_202011" # i.e. specify your build dir

setup_data() {
    # download data
    cd "pycode/memilio-epidata"
    python setup.py install
    python memilio/epidata/getPopulationData.py -o "$PWD/$data/pydata"
    python memilio/epidata/getDIVIData.py       -o "$PWD/$data/pydata" -s $startdate
    python memilio/epidata/getSimulationData.py -o "$PWD/$data/pydata" -s $startdate -m 7
    python memilio/epidata/getCaseData.py       -o "$PWD/$data/pydata" -s $startdate -m 7

    # run transformMobilityData.py with custom main (set directory to $PWD/$data/pydata, instead of using dict)
    python <<< "\
from memilio.epidata.transformMobilityData import *
directory = '$PWD/$data/mobility/'
updateMobility2022(directory, mobility_file='twitter_scaled_1252')
updateMobility2022(directory, mobility_file='commuter_migration_scaled')
# create federal states mobility matrix (not used in simulation for now)
createFederalStatesMobility(directory, mobility_file='twitter_scaled_1252')
createFederalStatesMobility(
    directory, mobility_file='commuter_migration_scaled')"

    cd $PWD
}

run_save() {
    n=$1
    save=$2
    result=$3

    mkdir -p $save
    mkdir -p $result
    
    printf "\nRunning: mpiexec -n $n $binary $data $save $result\n"
    mpiexec -n $n $binary $data $save $result
}

run_load() {
    n=$1
    save=$2
    result=$3

    mkdir -p $result
    
    printf "\nRunning: mpiexec -n $n $binary $save $result\n"
    mpiexec -n $n $binary $save $result
}

compare_numbers() {
    diff -qs $1 $2 >> $3
    if [ $? != 0 ]
    then
        bash $PWD/scripts/compare_numbers.sh $1 $2 $tolerance >> $3 
    fi
}

compare_h5() {
    # get text file from h5 and omit first line containing the files own path
    h5dump -w 0 -m "%.14g" $1 | tail -n +2 > "$results/1"
    h5dump -w 0 -m "%.14g" $2 | tail -n +2 > "$results/2"
    diff -qs "$results/1" "$results/2" >> $3
    different=$?
    echo "  with 1: $1" >> $3
    echo "       2: $2" >> $3
    if [ $different != 0 ]
    then
        compare_numbers "$results/1" "$results/2" $3
    fi
}

compare_recursive() {
    # create output directory
    mkdir -p "$3"
    # assert that dirs $1 and $2 have the same non empty file structure
    if ls $1 1> /dev/null && ls $2 1> /dev/null
        then
        ls1="$(ls $1)"
        if [[ $ls1 != "$(ls $2)" ]]
            then echo "Cannot compare directories with different file structures \"$1\" and \"$2\"." && return
        elif [[ $ls1 == "" ]]
            then echo "Cannot compare empty directories." && return
        fi
    fi

    # (over)write empty output file
    printf "" > "$3/diff.txt"
    # compare each file
    
    for path in $1/*
    do
        [[ -d $path ]] && continue # skip directories
        file=$(basename $path)
        if echo "$file" | grep -i ".h5" > /dev/null
            then compare_h5      "$1/$file" "$2/$file" "$3/diff.txt"
            else compare_numbers "$1/$file" "$2/$file" "$3/diff.txt"
        fi
    done

    for path in $1/*
    do
        ! [[ -d $path ]] && continue # skip files
        subdir=$(basename $path)
        compare_recursive "$1/$subdir" "$2/$subdir" "$3/$subdir"
    done
}

main() {
    # go to git root, just to be sure
    if git rev-parse --show-toplevel 1> /dev/null
    then cd $(git rev-parse --show-toplevel)
    else exit
    fi
    PWD=$(pwd)

    if ! [[ -f $virtual_env ]]
    then echo "Set virtual_env to an existing python venv." && exit
    else source $virtual_env
    fi
    # download and prepare data
    if [[ -d "$data/pydata/Germany" ]]
    then echo "Found Simulation data folder $data/pydata/Germany, skipping data setup"
    else setup_data
    fi

    # run simulations
    mkdir -p $saves
    # sequential
    run_save 1 "$saves/sequential" "$results/sequential/save"
    run_load 1 "$saves/sequential" "$results/sequential/load"
    # parallel
    run_save 5 "$saves/parallel" "$results/parallel/save"
    run_load 5 "$saves/parallel" "$results/parallel/load"

    #compare results
    compare_recursive "$results/sequential" "$results/parallel" "$output/diff/result"
    compare_recursive "$saves/sequential" "$saves/parallel" "$output/diff/save"
    # clean temp files
    rm "$results/1" "$results/2"

    # concat diff files into one
    cat $(find $output/diff/ -name diff.txt) > "$output/diff/gathered_diff.txt"
}

mkdir -p $output

if [ ! $? ]
then echo "Exiting" && exit
fi

if [[ -n "$(ls -A $output)" ]]
then echo "Use an empty \"output\" directory, as any preexisting files can lead to wrong results!"
fi

echo "Writing log to $output/log.txt"
printf "" > "$output/log.txt"
main &>> "$output/log.txt"
