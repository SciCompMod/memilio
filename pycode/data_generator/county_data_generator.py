import sys
import csv
import h5py
import random
import time
import epidemiology.secir as secir
import numpy as np
from epidemiology.secir import InfectionState as State


def get_county_model(county_data, day, dampings):
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

        model.populations[group, secir.Index_InfectionState(State.Exposed)] = county_data[i][day][1]
        model.populations[group, secir.Index_InfectionState(State.Carrier)] = county_data[i][day][2]
        model.populations[group, secir.Index_InfectionState(State.Infected)] = county_data[i][day][3]
        model.populations[group, secir.Index_InfectionState(State.Hospitalized)] = county_data[i][day][4]
        model.populations[group, secir.Index_InfectionState(State.ICU)] = county_data[i][day][5]
        model.populations[group, secir.Index_InfectionState(State.Recovered)] = county_data[i][day][6]
        model.populations[group, secir.Index_InfectionState(State.Dead)] = county_data[i][day][7]
        model.populations.set_difference_from_group_total_AgeGroup((group, secir.Index_InfectionState(State.Susceptible)), county_data[i][day][0])

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
    raw_baseline_matrices.append(open('../matrices/baseline_home.txt','r').read())
    raw_baseline_matrices.append(open('../matrices/baseline_work.txt','r').read())
    raw_baseline_matrices.append(open('../matrices/baseline_other.txt','r').read())
    raw_baseline_matrices.append(open('../matrices/baseline_school_pf_eig.txt','r').read())

    baseline_matrices = []
    for item in raw_baseline_matrices:
        baseline_matrices.append([[float(tmp) for tmp in row.split(' ')] for row in item.split('\n')])

    # sum over all matrices
    #baseline = [[baseline_matrices[0][ix][iy] + baseline_matrices[1][ix][iy] + baseline_matrices[2][ix][iy] + baseline_matrices[3][ix][iy] for iy in range(len(baseline_matrices[0][0]))] for ix in range(len(baseline_matrices[0]))]
    return np.array(baseline_matrices)


def getMinimumMatrix():
    raw_minimum_matrices = []
    raw_minimum_matrices.append(open('../matrices/minimum_home.txt','r').read())
    raw_minimum_matrices.append(open('../matrices/minimum_work.txt','r').read())
    raw_minimum_matrices.append(open('../matrices/minimum_other.txt','r').read())
    raw_minimum_matrices.append(open('../matrices/minimum_school_pf_eig.txt','r').read())

    minimum_matrices = []
    for item in raw_minimum_matrices:
        minimum_matrices.append([[float(tmp) for tmp in row.split(' ')] for row in item.split('\n')])

    #sum over all matrices
    #minimum = [[minimum_matrices[0][ix][iy] + minimum_matrices[1][ix][iy] + minimum_matrices[2][ix][iy] + minimum_matrices[3][ix][iy] for iy in range(len(minimum_matrices[0][0]))] for ix in range(len(minimum_matrices[0]))]
    return np.array(minimum_matrices)


def getRandomDampingMatrix():
    damping = np.zeros((6,6))
    for idx in range(6):
        for i in range(idx+1):
            damp = round(random.random(),4)
            damping[idx][i] = damp
            damping[i][idx] = damp

    return damping.tolist()


def getDampings(num_dampings, num_sim_days):
    dampings = []
    for num in range(num_dampings):
        #get random damping matrix
        dampings.append(getRandomDampingMatrix())
        #get random day for damping
        tmp = random.randint(0, num_sim_days-1)
        #check for doubles
        while(tmp in dampings):
            tmp = random.randint(0, num_sim_days-1)
        dampings.append(tmp)
    return dampings


def data_formatting(data, random_day, dampings, num_sim_days, results):
    #this list contains all input variables for the machine learning model

    all_county_samples = []
    for index in range(len(data)):
        x_input = list(dampings)
        for age_group in data[index]:
            rounded_compartment_values = [int(item) for item in age_group[:][random_day]]
            x_input = x_input + rounded_compartment_values
    # x_input format: damping, damping_day, damping, damping_day, ...... , [county][group]=> 8 values of compartments

        #this list contains all values from all compartments, age_groups, county's, on every day
        y_output = []
        #go trough all result runs (normaly its just one)
        for item in results:
            #iterate trough time steps
            for full_time_step in range(item.get_node(index).property.result.get_num_time_points()):
                #check whether the selected time point has no after point values. (No time points like: 1.2, 2.5, 3.4, just time points like 1.0, 2.0, 3.0)
                if not (item.get_node(index).property.result.get_time(full_time_step) - int(item.get_node(index).property.result.get_time(full_time_step))):
                    rounded_compartment_values = [int(item) for item in list(item.get_node(index).property.result.get_value(full_time_step))]
                    y_output = y_output + rounded_compartment_values
                    ##iterate trough all countys
                    #for county in range(len(data.keys())):
                    #    rounded_compartment_values = [int(item) for item in list(item.get_node(county).property.result.get_value(full_time_step))]
                    #    y_output = y_output + rounded_compartment_values

        all_county_samples.append(x_input + y_output)
    return all_county_samples




def single_simulation(data, day, dampings, num_sim_days):
    #basis graph, in the end this should simulate germany as country 
    graph = secir.SecirModelGraph()

    #add every county as a node in graph
    for index in range(len(data)):
      graph.add_node(index, get_county_model(data[index], day, dampings))

    #notes how to add edge, 
    #graph.add_edge(0, 1, 0.01 * np.ones(8*num_groups))
    #graph.add_edge(1, 0, 0.01 * np.ones(8*num_groups))



    study = secir.ParameterStudy(graph, t0 = 1, tmax = num_sim_days, dt = 0.5, num_runs = 1)
    return study.run()




def multiple_simulation_organizer(dataset_file_name, num_runs, num_dampings, num_sim_days):
    #load data
    data = h5py.File("Results_rki.h5","r")
    #print(list(data["1001"]["Group1"]))

    print("Start to format data")
    new_data = []
    group_names = ["Group1", "Group2", "Group3", "Group4", "Group5", "Group6"]
    for county in data.keys():
        new_data.append([data[county][name] for name in group_names])

    #change lines ~155 in code for data generation with multiple county output
    data = new_data
    #data = [list(np.sum(np.array(new_data), 0))]
    print("End formatting data")


    secure_save_after_num_runs = 5
    if(secure_save_after_num_runs > num_runs):
        secure_save_after_num_runs = num_runs

    output = []
    for save_circle in range(int(num_runs/secure_save_after_num_runs)):
        start = time.time()
        for run in range(secure_save_after_num_runs):

            #get random day from data
            #in results_rki.h5 are 180 samples, thus 0-179
            random_day = random.randint(0,179)

            #make dampings
            dampings = getDampings(num_dampings, num_sim_days)

            #run simulations
            results = single_simulation(data, random_day, dampings, num_sim_days)

            #format data to sample for machine learing
            formated_samples = data_formatting(data, random_day, dampings, num_sim_days, results)
            #add sample to output list
            output = output + formated_samples

        save(dataset_file_name, output, "a")
        output = []
        end = time.time()
        print("Secure save time point after: " + str(end-start) + " seconds")
        print(str(save_circle+1) + " of " + str(int(num_runs/secure_save_after_num_runs)) + "save time point circles")


def save(file_name, output, option):
    #store final data in csv file
    file = open(file_name, option)
    with file:
        writer = csv.writer(file)
        for row in output:
            writer.writerow(row)




if __name__ == "__main__":
    start = time.time()
    num_runs = int(sys.argv[1])
    dataset_file_name = sys.argv[2]
    num_dampings = int(sys.argv[3])
    num_sim_days = int(sys.argv[4])
    multiple_simulation_organizer(dataset_file_name=dataset_file_name, num_runs=num_runs , num_dampings=num_dampings, num_sim_days=num_sim_days)
    end = time.time()
    print(end - start)
