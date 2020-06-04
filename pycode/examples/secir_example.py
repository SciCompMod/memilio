from epidemiology.secir import ContactFrequencyMatrix, Damping, SecirParams, print_secir_params, simulate, StageTimes, Probabilities, Populations
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np


def plot_secir():
    """
    Simulates the secir model using the c++ code
    and plots the results
    """


    # set propability of infection
    inf_prob = 13/217
    use_dampings = True

    Mio = 1000000

    # load contact matrix
    contact_polymod = pd.read_csv('../epidemiology/contact_data/images/Polymod.csv', index_col=0).values

    damp_work = pd.read_csv('../epidemiology/contact_data/damp_work.csv', index_col=0).values
    damp_home = pd.read_csv('../epidemiology/contact_data/damp_home.csv', index_col=0).values
    damp_school = pd.read_csv('../epidemiology/contact_data/damp_school.csv', index_col=0).values
    damp_other = pd.read_csv('../epidemiology/contact_data/damp_other.csv', index_col=0).values


    contact_polymod = 0.5*(contact_polymod + contact_polymod.T)
    damp_work = 0.5 * (damp_work + damp_work.T)
    damp_home = 0.5 * (damp_home + damp_home.T)
    damp_school = 0.5 * (damp_school + damp_school.T)
    damp_other = 0.5 * (damp_other + damp_other.T)

    dampings = pd.read_csv('../epidemiology/contact_data/dampings.csv', index_col=0, header=[0, 1]).values
    #print(dampings)
    dampings = np.swapaxes(dampings.reshape(8,120,8),1,2)
    for i in range(120):
        dampings[:,:,i] = 0.5*(dampings[:,:,i] + dampings[:,:,i].T)

    print(np.array_equal(dampings[:, :, -2], dampings[:,:,-3]))


    num_groups = contact_polymod.shape[-1]


    times = StageTimes()
    times.set_incubation(5.2)  # R_2^(-1)+R_3^(-1)
    times.set_infectious_mild(6.)  # 4-14  (=R4^(-1))
    times.set_serialinterval(4.2)   # 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
    times.set_hospitalized_to_home(12.)  # 7-16 (=R5^(-1))
    times.set_home_to_hospitalized(5.)  # 2.5-7 (=R6^(-1))
    times.set_hospitalized_to_icu(2.)  # 1-3.5 (=R7^(-1))
    times.set_icu_to_home(8.)  # 5-16 (=R8^(-1))
    times.set_infectious_asymp(6.2)  # (=R9^(-1)=R_3^(-1)+0.5*R_4^(-1))
    times.set_icu_to_death(5.)  # 3.5-7 (=R5^(-1))



    sus_ORs=np.array([0.34, 0.67, 1.00, 1.00, 1.00, 1.00, 1.24, 1.47, 1.47]) # Odds ratios for relative susceptibility -- from https://science.sciencemag.org/content/early/2020/05/04/science.abb8001; 10-20 and 60-70 bins are the average across the ORs
    trans_ORs=np.array([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]) # Odds ratios for relative transmissibility -- no evidence of differences
    comorbidities=np.array([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00]) # Comorbidities by age -- set to 1 by default since already included in disease progression rates
    symp_probs=np.array([0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90])  # Overall probability of developing symptoms (based on https://www.medrxiv.org/content/10.1101/2020.03.24.20043018v1.full.pdf, scaled for overall symptomaticity)
    severe_probs=np.array([0.00050, 0.00165, 0.00720, 0.02080, 0.03430, 0.07650, 0.13280, 0.20655, 0.24570]) # Overall probability of developing severe symptoms (derived from Table 1 of https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-03-16-COVID19-Report-9.pdf)
    crit_probs=np.array([0.00003, 0.00008, 0.00036, 0.00104, 0.00216, 0.00933, 0.03639, 0.08923, 0.17420]) # Overall probability of developing critical symptoms (derived from Table 1 of https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-03-16-COVID19-Report-9.pdf)
    death_probs=np.array([0.00002, 0.00006, 0.00030, 0.00080, 0.00150, 0.00600, 0.02200, 0.05100, 0.09300]) # Overall probability of dying (https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf)

    '''probs = Probabilities()
    probs.set_asymp_per_infectious(0.09)  # 0.01-0.16
    probs.set_risk_from_symptomatic(0.25)  # 0.05-0.5
    probs.set_hospitalized_per_infectious(0.2)  # 0.1-0.35
    probs.set_icu_per_hospitalized(0.25)  # 0.15-0.4
    probs.set_dead_per_icu(0.3)  # 0.15-0.77'''
    probs = []
    people = []
    for i in range(num_groups):
        probs.append(Probabilities())
        probs[i].set_asymp_per_infectious(1-symp_probs[i])
        probs[i].set_risk_from_symptomatic(0.25)
        probs[i].set_hospitalized_per_infectious(severe_probs[i])
        probs[i].set_icu_per_hospitalized(crit_probs[i]/severe_probs[i])
        probs[i].set_dead_per_icu(death_probs[i]/(crit_probs[i]/severe_probs[i]))

        people.append(Populations())
        people[i].set_exposed_t0(100/num_groups)
        people[i].set_carrier_t0(50/num_groups)
        people[i].set_infectious_t0(50/num_groups)
        people[i].set_hospital_t0(10/num_groups)
        people[i].set_icu_t0(10/num_groups)
        people[i].set_recovered_t0(10/num_groups)
        people[i].set_dead_t0(0)


    # population numbers were taken from
    # https://de.statista.com/statistik/daten/studie/1365/umfrage/bevoelkerung-deutschlands-nach-altersgruppen/
    people[0].set_total_t0(4.66 * Mio)
    people[1].set_total_t0(8.93 * Mio)
    people[2].set_total_t0(11.31 * Mio)
    people[3].set_total_t0(10.84 * Mio)
    people[4].set_total_t0(13.9 * Mio)
    people[5].set_total_t0(10 * Mio)
    people[6].set_total_t0(9.5 * Mio)
    people[7].set_total_t0(13.88 * Mio)


    # set the params required for the simulation
    params = []
    for i in range(num_groups):
        params.append(SecirParams())
        params[i].times = times
        params[i].probabilities = probs[i]
        params[i].populations = people[i]
    print_secir_params(params)

    cont_freq_matrix = ContactFrequencyMatrix(num_groups)

    
    #  set contact rates and emulate some mitigations
    for i in range(num_groups):
        for j in range(0,num_groups):
            cont_freq_matrix.set_cont_freq(contact_polymod[i,j]*inf_prob, i, j)
            if use_dampings:
                for d in range(dampings.shape[2]):
                    cont_freq_matrix.add_damping(Damping(d, dampings[i,j,d]), i, j)
                #cont_freq_matrix.add_damping(Damping(25., damp_work[i,j]), i, j)
                #cont_freq_matrix.add_damping(Damping(35., damp_home[i,j]), i, j)
                #cont_freq_matrix.add_damping(Damping(40., damp_school[i,j]), i, j)
                #cont_freq_matrix.add_damping(Damping(50., damp_other[i,j]), i, j)

                #cont_freq_matrix.add_damping(Damping(10., 0.8), i, j)
                #cont_freq_matrix.add_damping(Damping(20., 0.9), i, j)
                #cont_freq_matrix.add_damping(Damping(30., 0.45), i, j)
                #cont_freq_matrix.add_damping(Damping(40., 0.25), i, j)
                #cont_freq_matrix.add_damping(Damping(70., 0.35), i, j)
                #cont_freq_matrix.add_damping(Damping(80., 0.45), i, j)
                #cont_freq_matrix.add_damping(Damping(100., 0.7), 0, 0)

    # run the simulation
    result = simulate(t0=0., tmax=120., dt=0.1, cont_freq_matrix=cont_freq_matrix, params=params)

    # sum over all groups
    group_data = np.zeros((len(result[0].t),8*num_groups))
    data = np.zeros((len(result[0].t),num_groups))
    for i in range(num_groups):
        group_data[:,0 + i*8] = result[i].nb_sus
        group_data[:,1 + i*8] = result[i].nb_exp
        group_data[:,2 + i*8] = result[i].nb_car
        group_data[:,3 + i*8] = result[i].nb_inf
        group_data[:,4 + i*8] = result[i].nb_hosp
        group_data[:,5 + i*8] = result[i].nb_icu
        group_data[:,6 + i*8] = result[i].nb_rec
        group_data[:,7 + i*8] = result[i].nb_dead

        data[:, 0] += result[i].nb_sus
        data[:, 1] += result[i].nb_exp
        data[:, 2] += result[i].nb_car
        data[:, 3] += result[i].nb_inf
        data[:, 4] += result[i].nb_hosp
        data[:, 5] += result[i].nb_icu
        data[:, 6] += result[i].nb_rec
        data[:, 7] += result[i].nb_dead



    fig, ax = plt.subplots()
    #ax.plot(result[0].t, data[:,0], label='#Suscepted')
    ax.plot(result[0].t, data[:,1], label='#Exposed')
    ax.plot(result[0].t, data[:,2], label='#Carrying')
    ax.plot(result[0].t, data[:,3], label='#Infected')
    ax.plot(result[0].t, data[:,4], label='#Hospitalzed')
    ax.plot(result[0].t, data[:,5], label='#In Intensive Care Units')
    #ax.plot(result[0].t, data[:,6], label='#Recovered')
    ax.plot(result[0].t, data[:,7], label='#Died')
    ax.set_title("SECIR model simulation")
    ax.set_xlabel("Time [days]")
    ax.legend()
    fig.savefig('Secir_Comix_Dampings_all.pdf')

    if True:
        compartiments = ['Suscepted', 'Exposed', 'Carrying', 'Infected', 'Hospitalized', 'ICU', 'Recoverd', 'Dead', 'time']
        groups = ['0-4', '5-17', '18-29', '30-39', '40-49', '50-59', '60-69', '70+']
        fig, ax = plt.subplots(2, 4, figsize=(30, 10))
        for i, title in zip(range(8), compartiments):
            for j, group in zip(range(num_groups), groups):
                ax[int(i % 2), int(np.floor(i / 2))].plot(result[j].t, group_data[:,j*num_groups+i], label=group)
            ax[int(i % 2), int(np.floor(i / 2))].set_title(title)
            ax[int(i % 2), int(np.floor(i / 2))].legend()
        fig.tight_layout()
        fig.savefig('Secir_Comix_Dampings_Groups.pdf')
    plt.show()
    plt.close()

    return data


if __name__ == "__main__":
    r1 = plot_secir()
    plt.show()