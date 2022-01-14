import epidemiology.secir as secir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime
import csv
import json
import time
import random
import sys, os
from _thread import start_new_thread

def load():
    # Define Comartment names
    compartments = ['Susceptible', 'Exposed', 'Carrying', 'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead', 'time']

    # define age groups
    groups = ['0-4', '5-14', '15-34', '35-59', '60-79', '80+']

    # processing population data to get the list with different count of each age group
    print("Load population data")
    populations = getPopulation()
    print("Population data successfully loaded")

    days = 100 # number of days to simulate
    t0 = 0
    dt = 0.1
    num_groups = len(groups)
    num_compartments = 8

    # normal contact matrix
    baseline_matrix = getBaselineMatrix()
    # minimum possible contact matrix
    minimum_matrix = getMinimumMatrix()


    # array of categories linked to input infos
    output_list = []

    print("Load and process RKI data")
    # load rki data
    with open('./data/all_county_age_rki_interpolated_ma.json') as f:
        rkidata = json.load(f)
    print("Loading successful")


    rki_deaths_dict = {}
    rki_recovered_dict = {}
    rki_infected_dict = {}
    for item in rkidata:
        rki_deaths_dict[item["Date"]] = [0 for i in range(num_groups)]
        rki_recovered_dict[item["Date"]] = [0 for i in range(num_groups)]
        rki_infected_dict[item["Date"]] = [0 for i in range(num_groups)]

    for item in rkidata:
        age_index = -1
        if(item["Age_RKI"] == "A00-A04"):
            age_index = 0
        elif(item["Age_RKI"] == "A05-A14"):
            age_index = 1
        elif(item["Age_RKI"] == "A15-A34"):
            age_index = 2
        elif(item["Age_RKI"] == "A35-A59"):
            age_index = 3
        elif(item["Age_RKI"] == "A60-A79"):
            age_index = 4
        elif(item["Age_RKI"] == "A80+"):
            age_index = 5

        rki_deaths_dict[item["Date"]][age_index] += item["Deaths"]
        rki_recovered_dict[item["Date"]][age_index] += item["Recovered"]
        rki_infected_dict[item["Date"]][age_index] += item["Confirmed"]

    rki_deaths = list(rki_deaths_dict.items())
    rki_recovered = list(rki_recovered_dict.items())

    rki_infected_tmp = list(rki_infected_dict.items())
    rki_infected = []
    for idx in range(len(rki_infected_tmp)):
        if(idx < 14):
            rki_infected.append((rki_infected_tmp[idx][0], [int(item) for item in rki_infected_tmp[idx][1]]))
        else:
            delta_day_before = [rki_infected_tmp[idx][1][i] - rki_infected_tmp[idx-1][1][i] for i in range(num_groups)]
            delta_14_days_before = [rki_infected_tmp[idx-13][1][i] - rki_infected_tmp[idx-14][1][i] for i in range(num_groups)]
            delta = [int(rki_infected[idx-1][1][i] + delta_day_before[i] - delta_14_days_before[i]) for i in range(num_groups)]
            rki_infected.append((rki_infected_tmp[idx][0], delta))

    #plt.plot([i for i in range(len(rki_infected))], [sum(item[1]) for item in rki_infected],'b', [i for i in range(len(rki_infected_tmp))], [sum(item[1]) for item in rki_infected_tmp], 'r')
    #plt.show()

    print("RKI data successfully preprocessed")


    print("Load DIVI data")
    # load divi intensivregister data (ICU/H data)
    with open('../../data/pydata/Germany/germany_divi.json') as f:
        divi_data_loaded = json.load(f)
    print("DIVI data successfully loaded")

    divi_data = [list(item.values()) for item in divi_data_loaded]

    # get first date from record of ICU data
    divi_data_start_date = divi_data[0][0]

    # get index from rki data which refers to the first divi data date (sychronization of datasets)
    rki_sychronized_index = [item[0] for item in rki_infected].index(divi_data_start_date)

    return groups, t0, days, dt, num_groups, populations, num_compartments, baseline_matrix, minimum_matrix, rki_sychronized_index, rki_infected, rki_recovered, rki_deaths, divi_data




def run(groups, t0, days, dt, num_groups, populations, num_compartments, baseline_matrix, minimum_matrix, rki_sychronized_index, rki_infected, rki_recovered, rki_deaths, divi_data, length, run_number, n_dampings):
    # final output list
    output_list = []
    for x in range(length):
        # list of the dampings
        # format dampings[n][0] => damping matrix, damping[n][1] => damping day
        dampings = []
        for num in range(n_dampings):
            damp = []
            # variation of damping
            damp.append(getRandomDampingMatrix(len(groups)))
            # variation of dampingday
            tmp = random.randint(0,days)
            # check for doubles
            while(any(tmp in sublist for sublist in dampings)):
                tmp = random.randint(0,days)

            damp.append(tmp)

            dampings.append(damp)

        # variation of start values (compartments and age groups) [data drom rki and divi]
        idx = random.randint(10,len(rki_infected)-rki_sychronized_index-1)
        idxrki = idx + rki_sychronized_index

        #print("idx:" + str(idx))
        rkidata = [rki_infected[idxrki][1], rki_recovered[idxrki][1], rki_deaths[idxrki][1], divi_data[idx][1]-divi_data[idx][2], divi_data[idx][2]]
        # Run Simulation
        result = secir.simulate(t0, days, dt, getmodel(dampings, num_groups, populations, rkidata, baseline_matrix, minimum_matrix))

        list_dampings = []
        for item in dampings:
            list_dampings.append(item[0])
            list_dampings.append(item[1])
        x_output = list_dampings + rkidata
        # x_output: damping_matrix[6,6], dampingday, rki_infected_aged[6], rki_recovered_aged[6], rki_deaths_aged[6], divi_data_hospital_total, divi_data_icu_total
        # all 1D-List's are not saved as list ==> every entry is saved in final output list
        y_output = processing_result(result, num_compartments, populations, num_groups, days)
        # all 1D-List's are not saved as list ==> every entry is saved in final output list

        output_list.append(x_output + y_output)

        if(x!=0 and x%5000==0):
            #status =["Index: " + str(i)+"\n" + "Estimated rest time: " + str((length-1-i)*(datetime.now() - tr)) + "\nCalculating time: " + str(datetime.now()-tr)]
            #save("Status_File.txt", status,"w")
            save("data/epi_data_"+ str(sys.argv[2]) +".csv", output_list,"a")
            output_list = []
    save("data/epi_data_"+ str(sys.argv[2]) +".csv", output_list,"a")
    #os.system("rm Status_File.txt")


def save(file_name, output, option):
    # store final data in csv file
    file = open(file_name,option)
    with file:
        writer = csv.writer(file)

        for row in output:
            writer.writerow(row)



def getBaselineMatrix():
    raw_baseline_matrices = []
    raw_baseline_matrices.append(open('./matrices/baseline_home.txt','r').read())
    raw_baseline_matrices.append(open('./matrices/baseline_work.txt','r').read())
    raw_baseline_matrices.append(open('./matrices/baseline_other.txt','r').read())
    raw_baseline_matrices.append(open('./matrices/baseline_school_pf_eig.txt','r').read())

    baseline_matrices = []
    for item in raw_baseline_matrices:
        baseline_matrices.append([[float(tmp) for tmp in row.split(' ')] for row in item.split('\n')])


    baseline = [[baseline_matrices[0][ix][iy] + baseline_matrices[1][ix][iy] + baseline_matrices[2][ix][iy] + baseline_matrices[3][ix][iy] for iy in range(len(baseline_matrices[0][0]))] for ix in range(len(baseline_matrices[0]))]
    return np.array(baseline)


def getMinimumMatrix():
    raw_minimum_matrices = []
    raw_minimum_matrices.append(open('./matrices/minimum_home.txt','r').read())
    raw_minimum_matrices.append(open('./matrices/minimum_work.txt','r').read())
    raw_minimum_matrices.append(open('./matrices/minimum_other.txt','r').read())
    raw_minimum_matrices.append(open('./matrices/minimum_school_pf_eig.txt','r').read())

    minimum_matrices = []
    for item in raw_minimum_matrices:
        minimum_matrices.append([[float(tmp) for tmp in row.split(' ')] for row in item.split('\n')])



    minimum = [[minimum_matrices[0][ix][iy] + minimum_matrices[1][ix][iy] + minimum_matrices[2][ix][iy] + minimum_matrices[3][ix][iy] for iy in range(len(minimum_matrices[0][0]))] for ix in range(len(minimum_matrices[0]))]
    return np.array(minimum)


def getRandomDampingMatrix(group_count):
    damping = np.zeros((group_count, group_count))
    for idx in range(group_count):
        for i in range(idx+1):
            damp = round(random.random(),4)
            damping[idx][i] = damp
            damping[i][idx] = damp

    return damping.tolist()



# processing popultation data ==> classification into different age groups
def getPopulation():
    # load population data from county_current_population.json
    with open('./data/adapted_county_current_population.json') as f:
        data = json.load(f)
    g0 = sum([item["0-4"] for item in data])
    g1 = sum([item["5-14"] for item in data])
    g2 = sum([item["15-34"] for item in data])
    g3 = sum([item["35-59"] for item in data])
    g4 = sum([item["60-79"] for item in data])
    g5 = sum([item["80+"] for item in data])
    return[g0,g1,g2,g3,g4,g5]






# classifying results
def processing_result(result, num_compartments, populations, num_groups, days):
    num_time_points = result.get_num_time_points()
    result_array = result.as_ndarray()
    t = result_array[0, :]
    group_data = np.transpose(result_array[1:, :])

    #age_data = np.zeros((num_time_points, num_compartments))
    age_data = []
    # age_data[age_group][day][compartment]
    for i in range(num_groups):
        #age_data += group_data[:,i*num_compartments:(i+1)*num_compartments]
        age_data.append(group_data[:,i*num_compartments:(i+1)*num_compartments])

    output = []

    for x in range(num_groups):
        for y in range(0,len(t),5):
            for z in range(num_compartments):
                output.append(round(age_data[x][y][z],1))

    return output



# parameter creation
# rkidata:  infected, recovered, deaths, icu
def getmodel(dampings, num_groups, populations, rkidata, baseline_matrix, minimum_matrix):
    # Initialize Parameters
    model = secir.SecirModel6()


    # Set parameters
    for i in range(num_groups):
        # Compartment transition duration
        model.parameters.times[i].set_incubation(5.2)
        model.parameters.times[i].set_infectious_mild(6)
        model.parameters.times[i].set_serialinterval(4.2)
        model.parameters.times[i].set_hospitalized_to_home(12)
        model.parameters.times[i].set_home_to_hospitalized(5)
        model.parameters.times[i].set_hospitalized_to_icu(2)
        model.parameters.times[i].set_icu_to_home(8)
        model.parameters.times[i].set_icu_to_death(5)

        # Initial number of peaople in each compartment
        model.populations.set(rkidata[0][i]/2, secir.AgeGroup6(i), secir.InfectionType.E)
        model.populations.set(rkidata[0][i]/2, secir.AgeGroup6(i), secir.InfectionType.C)
        model.populations.set(rkidata[0][i]/2, secir.AgeGroup6(i), secir.InfectionType.I)
        model.populations.set(int(rkidata[3]/15*(i+1)), secir.AgeGroup6(i), secir.InfectionType.H)
        model.populations.set(int(rkidata[4]/15*(i+1)), secir.AgeGroup6(i), secir.InfectionType.U)
        model.populations.set(rkidata[1][i], secir.AgeGroup6(i), secir.InfectionType.R)
        model.populations.set(rkidata[2][i], secir.AgeGroup6(i), secir.InfectionType.D)
        model.populations.set_difference_from_group_total_AgeGroup(populations[i], secir.AgeGroup6(i), secir.AgeGroup6(i), secir.InfectionType.S)

        # Compartment transition propabilities
        model.parameters.probabilities[i].set_infection_from_contact(1.0)
        model.parameters.probabilities[i].set_carrier_infectability(0.67)
        model.parameters.probabilities[i].set_asymp_per_infectious(0.09)
        model.parameters.probabilities[i].set_risk_from_symptomatic(0.25)
        model.parameters.probabilities[i].set_hospitalized_per_infectious(0.2)
        model.parameters.probabilities[i].set_icu_per_hospitalized(0.25)
        model.parameters.probabilities[i].set_dead_per_icu(0.3)


    # set contact frequency matrix
    model.parameters.get_contact_patterns().cont_freq_mat[0].baseline = baseline_matrix * 0.05
    model.parameters.get_contact_patterns().cont_freq_mat[0].minimum = minimum_matrix * 0.05

    # Define Damping on Contacts
    for item in dampings:
        model.parameters.get_contact_patterns().cont_freq_mat.add_damping(secir.Damping(item[0],item[1], 0, 0))

    # Apply mathematical constraints to parameters
    model.apply_constraints()

    return model



if __name__ == "__main__":
    groups, t0, days, dt, num_groups, populations, num_compartments, baseline_matrix, minimum_matrix, rki_sychronized_index, rki_infected, rki_recovered, rki_deaths, divi_data = load()

    length = int(sys.argv[1])
    run_number = int(sys.argv[2])
    n_dampings = int(sys.argv[3])
    run(groups, t0, days, dt, num_groups, populations, num_compartments, baseline_matrix, minimum_matrix, rki_sychronized_index, rki_infected, rki_recovered, rki_deaths, divi_data, length, run_number, n_dampings)
