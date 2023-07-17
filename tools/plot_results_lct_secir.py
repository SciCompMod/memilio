#############################################################################
# Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
#
# Authors: Lena Ploetzke
#
# Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#############################################################################

import h5py
import os
import pandas as pd
import matplotlib.pyplot as plt

import memilio.epidata.getDataIntoPandasDataFrame as gd



def plot_lct_result(file, compartment_idx=range(8)):
    # Define compartments
    secir_dict = {0:'Susceptible', 1:'Exposed',2:'Carrier', 3:'Infected', 4:'Hospitalized',
              5:'ICU', 6:'Recovered', 7:'Death'}

    input_file = os.path.join(os.getcwd(), str(file))
    h5file = h5py.File(input_file + '.h5', 'r')
    if(len(list(h5file.keys()))>1):
        raise gd.DataError("File should contain one dataset.")
    if(len(list(h5file[list(h5file.keys())[0]].keys()))>3):
        raise gd.DataError("Expected only one group.")


    data=h5file[list(h5file.keys())[0]]
    dates = data['Time'][:]
    total = data['Total'][:,compartment_idx]

    plt.figure('Vergleich f_{beta} der Simulationsmodelle')
    plt.plot(dates, total, linewidth=1.0)
    legendplot=[]
    for i in compartment_idx:
        legendplot.append(secir_dict[i])
    plt.legend(legendplot, fontsize=14)
    plt.xlabel('Zeit', fontsize=14)
    plt.ylabel('Anzahl Personen', fontsize=10)
    plt.grid(True, linestyle='--')
    plt.show()
    h5file.close()
    
    
    

if __name__ == '__main__':
    arr=list(range(1,5))
    arr.append(7)
    plot_lct_result("result_lct",arr)