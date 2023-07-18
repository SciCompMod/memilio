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

def get_subcompartments():
    # define the used subcompartments for simulation of the LCT model 
    vec_subcompartments=[1,2,3,1,1,5,1,1] 
    lct_secir_dict= {0:'S', 1:'E1', 2:'E2', 3:'C1', 4:'C2',
              5:'C3', 6:'I', 7:'H', 8:'U1', 9:'U2', 10:'U3', 11:'U4', 12:'U5', 13:'R', 14: 'D'}
    return (vec_subcompartments, lct_secir_dict)

def plot_lct_subcompartments(file, save = True):
    (vec_subcompartments, lct_secir_dict)=get_subcompartments()

    input_file = os.path.join(os.getcwd(), str(file))
    h5file = h5py.File(input_file + '.h5', 'r')

    if(len(list(h5file.keys()))>1):
        raise gd.DataError("File should contain one dataset.")
    if(len(list(h5file[list(h5file.keys())[0]].keys()))>3):
        raise gd.DataError("Expected only one group.")
    
    data=h5file[list(h5file.keys())[0]]
    dates = data['Time'][:]
    # As there should be only one Group, total is the simulation result 
    total = data['Total'][:,:]
    if(total.shape[1]!=sum(vec_subcompartments)):
        raise gd.DataError("Expected a different number of subcompartments.")
    
    fig, axs = plt.subplots(4,2, sharex='all', num='Subcompartments LCT SECIR model fictional scenario')
    for i in range(len(vec_subcompartments)):
        idx_start=sum(vec_subcompartments[0:i])
        axs[int(i/2),i%2].plot(dates,total[:,idx_start:idx_start+vec_subcompartments[i]])
        axs[int(i/2),i%2].grid(True, linestyle='--')
        legendplot=[]
        for j in range(idx_start,idx_start+vec_subcompartments[i]):
            legendplot.append(lct_secir_dict[j])
        axs[int(i/2),i%2].legend(legendplot, fontsize=10)
        axs[int(i/2),i%2].set_ylim(bottom=0)
        axs[int(i/2),i%2].set_xlim(left=0)
    fig.supxlabel('Time')
    fig.supylabel('Number of persons')
    
    if save:
        if not os.path.isdir('Plots'):
            os.makedirs('Plots')
        fig.savefig('Plots/lct_secir_fictional_subcompartments.png', bbox_inches='tight', dpi=500)
    plt.show()
    h5file.close()
    


def plot_lct_result(file, compartment_idx=range(8), save = True):
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
    # As there should be only one Group, total is the simulation result 
    total = data['Total'][:,compartment_idx]

    subcompartments=get_subcompartments()[0]
    heading_sub=', '.join(str(sub) for sub in subcompartments)
    plt.figure('SECIR LCT model, subcompartments:['+heading_sub+'], fictional scenario')
    plt.plot(dates, total, linewidth=1.0)

    legendplot=[]
    for i in compartment_idx:
        legendplot.append(secir_dict[i])
    plt.legend(legendplot, fontsize=14)

    plt.xlabel('Time', fontsize=14)
    plt.ylabel('Number of persons', fontsize=10)
    plt.ylim(bottom=0)
    plt.xlim(left=0)
    plt.grid(True, linestyle='--')
    if save:
        if not os.path.isdir('Plots'):
            os.makedirs('Plots')
        plt.savefig('Plots/lct_secir_fictional.png', bbox_inches='tight', dpi=500)
    plt.show()
    h5file.close()

def compare_results(files, legendplot ,save = True):
    secir_dict = {0:'Susceptible', 1:'Exposed',2:'Carrier', 3:'Infected', 4:'Hospitalized',
              5:'ICU', 6:'Recovered', 7:'Death'}
    
    fig, axs = plt.subplots(4,2, sharex='all', num='Compare files')

    for file in range(len(files)):
        input_file = os.path.join(os.getcwd(), str(files[file]))
        h5file = h5py.File(input_file + '.h5', 'r')

        if(len(list(h5file.keys()))>1):
            raise gd.DataError("File should contain one dataset.")
        if(len(list(h5file[list(h5file.keys())[0]].keys()))>3):
            raise gd.DataError("Expected only one group.")
    
        data=h5file[list(h5file.keys())[0]]
        dates = data['Time'][:]
        # As there should be only one Group, total is the simulation result 
        total = data['Total'][:,:]
        if(total.shape[1]!=8):
            raise gd.DataError("Expected a different number of subcompartments.")
    
        for i in range(8):
            axs[int(i/2),i%2].plot(dates,total[:,i],label=legendplot[file])
        
        h5file.close()

    for i in range(8):
        axs[int(i/2),i%2].set_title(secir_dict[i], fontsize=8)
        axs[int(i/2),i%2].set_ylim(bottom=0)
        axs[int(i/2),i%2].set_xlim(left=0)
        axs[int(i/2),i%2].grid(True, linestyle='--')
        axs[int(i/2),i%2].legend(fontsize=8)

    fig.supxlabel('Time')
    fig.supylabel('Number of persons')

    if save:
        if not os.path.isdir('Plots'):
            os.makedirs('Plots')
        fig.savefig('Plots/lct_secir_fictional_subcompartments.png', bbox_inches='tight', dpi=500)
    plt.show() 
    
if __name__ == '__main__':
    arr=list(range(1,5))
    arr.append(7)
    plot_lct_result("result_lct",arr)

    compare_results(["result_lct","result_ode"],legendplot=list(["LCT","ODE"]), save =True)
    