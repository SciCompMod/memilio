#!/bin/bash
# the regex does not accept commas, e.g. 1,000.0 or 0,5 are considered two numbers
number_regex="[+-]?[0-9]+([\.][0-9]+)?([eE][+-]?[0-9]+)?"
eps="1e-16"
if (($# < 2))
then
    printf "Compare numbers in order of appearance from two given files, up to some tolerance (by default 1e-16).\n"
    printf "Usage:\n$0 <file 1> <file 2> <optional: tolerance>\n"
    exit
fi
if (($# >= 3))
then eps="$3"
fi

numbers1=( $(grep -oE $number_regex $1) )
numbers2=( $(grep -oE $number_regex $2) )

python3 <<< "
import numpy as np

eps = float($eps)
eps = np.abs(eps)

numbers1 = np.array(str('${numbers1[*]}').split(), dtype=float)
numbers2 = np.array(str('${numbers2[*]}').split(), dtype=float)

if len(numbers1) != len(numbers2):
    print('  Files do not have the same amount of numbers!')
    print('  File 1: Found {} numbers.'.format(len(numbers1)))
    print('  File 2: Found {} numbers.'.format(len(numbers2)))
else:
    mask = np.where(numbers1 != numbers2)[0]
    if len(mask) == 0:
        print('  All {} numbers in Files \'$1\' and \'$2\' match exactly.'.format(len(numbers1)))
    else:
        comp = (np.abs(numbers1[mask] - numbers2[mask]) < eps)
        if np.all(comp):
            print('  All {} numbers in Files \'$1\' and \'$2\' match, with a tolerance of {}.'.format(len(numbers1), eps))
        else:
            print('  Mismatch between Files \'$1\' and \'$2\':')
            mask = mask[np.where(~comp)[0]]
            for m in mask:
                print('    {}: {} != {}'.format(m, numbers1[m], numbers2[m]))
"