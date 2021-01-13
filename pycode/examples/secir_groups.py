from epidemiology.secir import (UncertainContactMatrix, ContactMatrix, Damping, SecirModel8, simulate,
                                AgeGroup8, InfectionType)
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import sys
from datetime import datetime
import datetime as dt

from dampings import create_dampings


def plot_secir():
    """
    Simulates the secir model using the c++ code
    and plots the results
    """


    # set propability of infection
    # number of infected/number of close contacts -- from https://www.sciencedirect.com/science/article/pii/S1201971220302502#bib0045
    inf_prob =132/2147


    days = 120
    name = '_contact'
    use_dampings = True

    Mio = 1000000

    # load contact matrix
    contact_polymod = pd.read_csv('../../data/contact_data/images/Polymod.csv', index_col=0).values


    contact_polymod = 0.5*(contact_polymod + contact_polymod.T)
    
    num_groups = contact_polymod.shape[-1]
    assert num_groups == 8
    num_compartments = int(InfectionType.Count)

    create_dampings(path='../../data/contact_data/', days=days)
    dampings = pd.read_csv('../../data/contact_data/dampings.csv', index_col=0, header=[0, 1]).values
    dampings = np.swapaxes(dampings.reshape(num_groups,days,num_groups),1,2)
    for i in range(days):
        dampings[:,:,i] = 0.5*(dampings[:,:,i] + dampings[:,:,i].T)

    sus_ORs=np.array([0.34, 0.67, 1.00, 1.00, 1.00, 1.00, 1.24, 1.47]) # Odds ratios for relative susceptibility -- from https://science.sciencemag.org/content/early/2020/05/04/science.abb8001; 10-20 and 60-70 bins are the average across the ORs
    symp_probs=np.array([0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, (0.85 + 0.90)/2])  # Overall probability of developing symptoms (based on https://www.medrxiv.org/content/10.1101/2020.03.24.20043018v1.full.pdf, scaled for overall symptomaticity)
    severe_probs=np.array([0.00050, 0.00165, 0.00720, 0.02080, 0.03430, 0.07650, 0.13280, (0.20655 + 0.24570)/2]) # Overall probability of developing severe symptoms (derived from Table 1 of https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-03-16-COVID19-Report-9.pdf)
    crit_probs=np.array([0.00003, 0.00008, 0.00036, 0.00104, 0.00216, 0.00933, 0.03639, (0.08923 + 0.17420)/2]) # Overall probability of developing critical symptoms (derived from Table 1 of https://www.imperial.ac.uk/media/imperial-college/medicine/mrc-gida/2020-03-16-COVID19-Report-9.pdf)
    death_probs=np.array([0.00002, 0.00006, 0.00030, 0.00080, 0.00150, 0.00600, 0.02200, (0.05100 + 0.09300)/2]) # Overall probability of dying (https://www.imperial.ac.uk/media/imperial-college/medicine/sph/ide/gida-fellowships/Imperial-College-COVID19-NPI-modelling-16-03-2020.pdf)

    # population numbers were taken from
    # https://de.statista.com/statistik/daten/studie/1365/umfrage/bevoelkerung-deutschlands-nach-altersgruppen/
    # https://www.bpb.de/nachschlagen/zahlen-und-fakten/soziale-situation-in-deutschland/61538/altersgruppen
    num_pop = [3.88, 9.71, 12.53, 9.8, 13.7, 11.69, 9.0, 12.4]

    '''probs = Probabilities()
    probs.set_infection_from_contact(1.0)
    probs.set_asymp_per_infectious(0.09)  # 0.01-0.16
    probs.set_risk_from_symptomatic(0.25)  # 0.05-0.5
    probs.set_hospitalized_per_infectious(0.2)  # 0.1-0.35
    probs.set_icu_per_hospitalized(0.25)  # 0.15-0.4
    probs.set_dead_per_icu(0.3)  # 0.15-0.77'''

    # set the params required for the simulation
    model = SecirModel8()
    for i in range(num_groups):
        model.parameters.times[i].set_incubation(5.2)  # R_2^(-1)+R_3^(-1)
        model.parameters.times[i].set_infectious_mild(6.)  # 4-14  (=R4^(-1))
        model.parameters.times[i].set_serialinterval(4.2)  # 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        model.parameters.times[i].set_hospitalized_to_home(12.)  # 7-16 (=R5^(-1))
        model.parameters.times[i].set_home_to_hospitalized(5.)  # 2.5-7 (=R6^(-1))
        model.parameters.times[i].set_hospitalized_to_icu(2.)  # 1-3.5 (=R7^(-1))
        model.parameters.times[i].set_icu_to_home(8.)  # 5-16 (=R8^(-1))
        model.parameters.times[i].set_infectious_asymp(6.2)  # (=R9^(-1)=R_3^(-1)+0.5*R_4^(-1))
        model.parameters.times[i].set_icu_to_death(5.)  # 3.5-7 (=R5^(-1))

        model.parameters.probabilities[i].set_infection_from_contact(sus_ORs[i])
        model.parameters.probabilities[i].set_carrier_infectability(0.67)
        model.parameters.probabilities[i].set_asymp_per_infectious(1-symp_probs[i])
        model.parameters.probabilities[i].set_risk_from_symptomatic(0.25)
        model.parameters.probabilities[i].set_hospitalized_per_infectious(severe_probs[i])
        model.parameters.probabilities[i].set_icu_per_hospitalized(crit_probs[i]/severe_probs[i])
        model.parameters.probabilities[i].set_dead_per_icu(death_probs[i]/(crit_probs[i]/severe_probs[i]))

        model.populations.set(100*num_pop[i]/sum(num_pop), AgeGroup8(i), InfectionType.E)
        model.populations.set( 50*num_pop[i]/sum(num_pop), AgeGroup8(i), InfectionType.C)
        model.populations.set( 50*num_pop[i]/sum(num_pop), AgeGroup8(i), InfectionType.I)
        model.populations.set( 10*num_pop[i]/sum(num_pop), AgeGroup8(i), InfectionType.H)
        model.populations.set( 10*num_pop[i]/sum(num_pop), AgeGroup8(i), InfectionType.U)
        model.populations.set( 10*num_pop[i]/sum(num_pop), AgeGroup8(i), InfectionType.R)
        model.populations.set(  0*num_pop[i]/sum(num_pop), AgeGroup8(i), InfectionType.D)

    model.populations.set_total(Mio)
    
    #  set contact rates and emulate some mitigations
    for i in range(num_groups):
        for j in range(0,num_groups):
            model.parameters.get_contact_patterns().get_cont_freq_mat().set_cont_freq(contact_polymod[i,j]*inf_prob, i, j)
            # params.get_contact_patterns().get_cont_freq_mat().set_cont_freq(0.62*num_pop[i]/sum(num_pop), i, j)
            if use_dampings:
                for d in range(dampings.shape[2]):
                    model.parameters.get_contact_patterns().get_cont_freq_mat().add_damping(Damping(d, dampings[i,j,d]*1.06), i, j)
                #  params.get_contact_patterns().get_cont_freq_mat().add_damping(Damping(25., damp_work[i,j]), i, j)
                #  params.get_contact_patterns().get_cont_freq_mat().add_damping(Damping(35., damp_home[i,j]), i, j)
                #  params.get_contact_patterns().get_cont_freq_mat().add_damping(Damping(40., damp_school[i,j]), i, j)
                #  params.get_contact_patterns().get_cont_freq_mat().add_damping(Damping(50., damp_other[i,j]), i, j)

                # params.get_contact_patterns().get_cont_freq_mat().add_damping(Damping(25., 0.8), i, j)
                # params.get_contact_patterns().get_cont_freq_mat().add_damping(Damping(35., 0.9), i, j)
                # params.get_contact_patterns().get_cont_freq_mat().add_damping(Damping(40., 0.45), i, j)
                # params.get_contact_patterns().get_cont_freq_mat().add_damping(Damping(50., 0.25), i, j)
                # params.get_contact_patterns().get_cont_freq_mat().add_damping(Damping(70., 0.35), i, j)
                # params.get_contact_patterns().get_cont_freq_mat().add_damping(Damping(80., 0.45), i, j)
                # params.get_contact_patterns().get_cont_freq_mat().add_damping(Damping(100., 0.7), 0, 0)

    # run the simulation
    result = simulate(t0=0., tmax=days, dt=0.1, model=model)

    num_time_points = result.get_num_time_points()
    result_array = result.as_ndarray()
    t = result_array[0, :]
    group_data = np.transpose(result_array[1:, :])

    #sum over all groups
    data = np.zeros((num_time_points,num_compartments))
    for i in range(num_groups):
        data += group_data[:, i*num_compartments:(i + 1) * num_compartments]

    datelist = np.array(pd.date_range(datetime(2020, 1, 27), periods=days, freq='D').strftime('%m-%d').tolist())
    datelist2 = np.array(pd.date_range(datetime(2020, 2, 3), periods=days, freq='D').strftime('%m-%d').tolist())

    tick_range = (np.arange(int(days / 10) +1 ) * 10)
    tick_range[-1] -=1
    fig, ax = plt.subplots()
    #ax.plot(t, data[:,0], label='#Suscepted')
    ax.plot(t, data[:,1], label='#Exposed')
    ax.plot(t, data[:,2], label='#Carrying')
    ax.plot(t, data[:,3], label='#Infected')
    ax.plot(t, data[:,4], label='#Hospitalzed')
    ax.plot(t, data[:,5], label='#In Intensive Care Units')
    #ax.plot(t, data[:,6], label='#Recovered')
    ax.plot(t, data[:,7], label='#Died')
    ax.set_title("SECIR model simulation")
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range],rotation=45)
    ax.legend()
    fig.tight_layout
    fig.savefig('Secir_all' + name + '.pdf')

    ind_list = []
    for i in range(days):
        ind_list.append(np.argmin(np.abs(t-i)))

    new_sus = []
    for i in range(1,len(ind_list)):
        new_sus.append((data[ind_list[i-1],0] - data[ind_list[i], 0])/(t[ind_list[i]]-t[ind_list[i-1]]))

    fig, ax = plt.subplots()
    ax.plot(range(days-1), new_sus)
    ax.set_ylabel('Number of new cases')
    ax.set_title('Number of new Infections')
    ax.set_xticks(tick_range)
    ax.set_xticklabels(datelist[tick_range], rotation=45)
    fig.tight_layout
    fig.savefig('new_infections' + name + '.pdf')


    if True:
        compartments = ['Suscepted', 'Exposed', 'Carrying', 'Infected', 'Hospitalized', 'ICU', 'Recovered', 'Dead', 'time']
        groups = ['0-4', '5-17', '18-29', '30-39', '40-49', '50-59', '60-69', '70+']
        fig, ax = plt.subplots(4, 2, figsize=(12, 15))
        for i, title in zip(range(num_compartments), compartments):
            for j, group in zip(range(num_groups), groups):
                ax[int(np.floor(i / 2)),int(i % 2)].plot(t, group_data[:,j*num_groups+i], label=group)
            ax[int(np.floor(i / 2)),int(i % 2)].set_title(title)
            ax[int(np.floor(i / 2)), int(i % 2)].legend()

            ax[int(np.floor(i / 2)), int(i % 2)].set_xticks(tick_range)
            ax[int(np.floor(i / 2)), int(i % 2)].set_xticklabels(datelist[tick_range], rotation=45)
        plt.subplots_adjust(hspace=0.3, bottom=0.05, top=0.95)
        #fig.tight_layout()
        fig.savefig('Secir_Groups' + name + '.pdf')


        fig, ax = plt.subplots(4, 2, figsize=(12, 15))
        for i, title in zip(range(num_compartments), compartments):
            ax[int(np.floor(i / 2)), int(i % 2)].plot(t, data[:, i])
            ax[int(np.floor(i / 2)), int(i % 2)].set_title(title)

            ax[int(np.floor(i / 2)), int(i % 2)].set_xticks(tick_range)
            ax[int(np.floor(i / 2)), int(i % 2)].set_xticklabels(datelist[tick_range], rotation=45)
        ax[0,0].set_ylim(0, 90 * Mio)
        plt.subplots_adjust(hspace=0.3, bottom=0.05, top=0.95)
        # fig.tight_layout()
        fig.savefig('Secir_all_parts' + name + '.pdf')
    plt.show()
    plt.close()

    return data


if __name__ == "__main__":
    r1 = plot_secir()
    plt.show()
