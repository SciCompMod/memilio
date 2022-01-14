import epidemiology.secir as secir
from epidemiology.secir import InfectionState as State
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import csv
import h5py
import json
import time
import random
import sys, os
from _thread import start_new_thread



def single_simulation_run(day, dampings, num_sim_days):
    #load data
    data = h5py.File("data_generator/Results_rki.h5","r")
    #print(list(data["1001"]["Group1"]))

    print("Start to format data")
    new_data = []
    group_names = ["Group1", "Group2", "Group3", "Group4", "Group5", "Group6"]
    for county in data.keys():
        new_data.append([data[county][name][day] for name in group_names])

    #change lines ~155 in code for data generation with multiple county output
    data = [new_data[150]]
    #data = [list(np.sum(np.array(new_data), 0))]
    print("End formatting data")


    graph = secir.SecirModelGraph()

    #add every county as a node in graph
    for index in range(len(data)):
      graph.add_node(index, get_county_model(data[index], dampings))

    t1 = time.time()
    sim_result = secir.ParameterStudy(graph, t0 = 1, tmax = num_sim_days, dt = 0.5, num_runs = 1)
    sim_time = time.time() - t1
    results = data_formatting(data, dampings, num_sim_days, sim_result.run())
    return results, sim_time





def data_formatting(data, dampings, num_sim_days, results):
    #this list contains all input variables for the machine learning model
    x_input = list(dampings)

    for county in data:
        for age_group in county:
            rounded_compartment_values = [int(item) for item in age_group]
            x_input = x_input + rounded_compartment_values
    # x_input format: damping, damping_day, damping, damping_day, ...... , [county][group]=> 8 values of compartments

    #this list contains all values from all compartments, age_groups, county's, on every day
    y_output = []
    #go trough all result runs (normaly its just one)
    for item in results:
        #iterate trough time steps
        for full_time_step in range(item.get_node(0).property.result.get_num_time_points()):
            #check whether the selected time point has no after point values. (No time points like: 1.2, 2.5, 3.4, just time points like 1.0, 2.0, 3.0)
            if not (item.get_node(0).property.result.get_time(full_time_step) - int(item.get_node(0).property.result.get_time(full_time_step))):
                rounded_compartment_values = [int(item) for item in list(item.get_node(0).property.result.get_value(full_time_step))]
                y_output = y_output + rounded_compartment_values
                ##iterate trough all countys
                #for county in range(len(data.keys())):
                #    rounded_compartment_values = [int(item) for item in list(item.get_node(county).property.result.get_value(full_time_step))]
                #    y_output = y_output + rounded_compartment_values

    sample = [x_input, y_output]
    return sample




# parameter creation
# rkidata:  infected, recovered, deaths, icu
def get_county_model(county_data, dampings):
    #setup basic parameters
    num_groups = 6
    model = secir.SecirModel(num_groups)


    for i in range(num_groups):
        group = secir.AgeGroup(i)

        model.parameters.IncubationTime[group] = 5.2
        model.parameters.InfectiousTimeMild[group] = 7
        model.parameters.SerialInterval[group] = 4.25
        model.parameters.HospitalizedToHomeTime[group] = [5, 5, 6, 8, 10, 15][i]
        model.parameters.HomeToHospitalizedTime[group] = [10.5, 10.5, 10.5, 6, 6, 6][i]
        model.parameters.HospitalizedToICUTime[group] = 5
        model.parameters.ICUToHomeTime[group] = [7, 7, 7, 17.5, 17.5, 12.5][i]
        model.parameters.ICUToDeathTime[group] = [6, 6, 6, 16.5, 16.5, 11][i]

        model.populations[group, secir.Index_InfectionState(State.Exposed)] = county_data[i][1]
        model.populations[group, secir.Index_InfectionState(State.Carrier)] = county_data[i][2]
        model.populations[group, secir.Index_InfectionState(State.Infected)] = county_data[i][3]
        model.populations[group, secir.Index_InfectionState(State.Hospitalized)] = county_data[i][4]
        model.populations[group, secir.Index_InfectionState(State.ICU)] = county_data[i][5]
        model.populations[group, secir.Index_InfectionState(State.Recovered)] = county_data[i][6]
        model.populations[group, secir.Index_InfectionState(State.Dead)] = county_data[i][7]
        model.populations.set_difference_from_group_total_AgeGroup((group, secir.Index_InfectionState(State.Susceptible)), county_data[i][0])

        f = 1.4
        model.parameters.InfectionProbabilityFromContact[group].set_distribution(secir.ParameterDistributionUniform([0.02 * f, 0.05 * f,0.05 * f, 0.05 * f, 0.08 * f, 0.15 * f][i], [0.04 * f, 0.07 * f, 0.07 * f, 0.07 * f, 0.1 * f, 0.2 * f][i]))
        model.parameters.AsymptoticCasesPerInfectious[group] = [0.25, 0.25, 0.2, 0.2, 0.2, 0.2][i]
        model.parameters.RiskOfInfectionFromSympomatic[group] = 0.5
        model.parameters.HospitalizedCasesPerInfectious[group] = [0.0075, 0.0075, 0.019, 0.0625, 0.0165, 0.225][i]
        model.parameters.ICUCasesPerHospitalized[group] = [0.075, 0.075, 0.075, 0.15, 0.3, 0.4][i]
        model.parameters.DeathsPerHospitalized[group] = [0.05, 0.05, 0.14, 0.14, 0.4, 0.6][i]


    baseline = getBaselineMatrix()
    minimum =  getMinimumMatrix()

    model.parameters.ContactPatterns.cont_freq_mat = secir.ContactMatrixGroup(4,num_groups)
    model.parameters.ContactPatterns.cont_freq_mat[0].baseline = baseline[0]
    model.parameters.ContactPatterns.cont_freq_mat[0].minimum = minimum[0]
    model.parameters.ContactPatterns.cont_freq_mat[1].baseline = baseline[1]
    model.parameters.ContactPatterns.cont_freq_mat[1].minimum = minimum[1]
    model.parameters.ContactPatterns.cont_freq_mat[2].baseline = baseline[2]
    model.parameters.ContactPatterns.cont_freq_mat[2].minimum = minimum[2]
    model.parameters.ContactPatterns.cont_freq_mat[3].baseline = baseline[3]
    model.parameters.ContactPatterns.cont_freq_mat[3].minimum = minimum[3]

    #add dampings (list dampings: [damping_matrix, damping_day, damping_matrix, damping_day, .....])
    for i in range(int(len(dampings)/2)):
        model.parameters.ContactPatterns.cont_freq_mat.add_damping(secir.Damping(dampings[i*2], dampings[i*2+1]))

    #print(model.parameters.ContactPatterns.cont_freq_mat[1].baseline)
    #process the result of one run
    #parameter_study.c = 0


    #study the effect of different infection rates

    model.apply_constraints()

    return model


def getBaselineMatrix():
    raw_baseline_matrices = []
    raw_baseline_matrices.append(open('./matrices/baseline_home.txt','r').read())
    raw_baseline_matrices.append(open('./matrices/baseline_work.txt','r').read())
    raw_baseline_matrices.append(open('./matrices/baseline_other.txt','r').read())
    raw_baseline_matrices.append(open('./matrices/baseline_school_pf_eig.txt','r').read())

    baseline_matrices = []
    for item in raw_baseline_matrices:
        baseline_matrices.append([[float(tmp) for tmp in row.split(' ')] for row in item.split('\n')])

    # sum over all matrices
    #baseline = [[baseline_matrices[0][ix][iy] + baseline_matrices[1][ix][iy] + baseline_matrices[2][ix][iy] + baseline_matrices[3][ix][iy] for iy in range(len(baseline_matrices[0][0]))] for ix in range(len(baseline_matrices[0]))]
    return np.array(baseline_matrices)


def getMinimumMatrix():
    raw_minimum_matrices = []
    raw_minimum_matrices.append(open('./matrices/minimum_home.txt','r').read())
    raw_minimum_matrices.append(open('./matrices/minimum_work.txt','r').read())
    raw_minimum_matrices.append(open('./matrices/minimum_other.txt','r').read())
    raw_minimum_matrices.append(open('./matrices/minimum_school_pf_eig.txt','r').read())

    minimum_matrices = []
    for item in raw_minimum_matrices:
        minimum_matrices.append([[float(tmp) for tmp in row.split(' ')] for row in item.split('\n')])

    #sum over all matrices
    #minimum = [[minimum_matrices[0][ix][iy] + minimum_matrices[1][ix][iy] + minimum_matrices[2][ix][iy] + minimum_matrices[3][ix][iy] for iy in range(len(minimum_matrices[0][0]))] for ix in range(len(minimum_matrices[0]))]
    return np.array(minimum_matrices)
