from epidemiology.secir import Damping, SeirParam, print_seir_params, simulate, SecirResult
from matplotlib import pyplot as plt


def plot_secir():
    """
    Simulates the secir model using the c++ code
    and plots the results
    """

    tinc = 5.2 #  R_2^(-1)+R_3^(-1)
    tinfmild = 6 #  4-14  (=R4^(-1))
    tserint = 4.2 #  4-4.4 // R_2^(-1)+0.5*R_3^(-1)
    thosp2home = 12 #  7-16 (=R5^(-1))
    thome2hosp = 5 #  2.5-7 (=R6^(-1))
    thosp2icu = 2 #  1-3.5 (=R7^(-1))
    ticu2home = 8 #  5-16 (=R8^(-1))
    tinfasy = 6.2 #  (=R9^(-1)=R_3^(-1)+0.5*R_4^(-1))
    ticu2death = 5 # 3.5-7 (=R5^(-1))

    cont_freq = 0.5 #  0.2-0.75
    alpha = 0.09 #  0.01-0.16
    beta = 0.25 #  0.05-0.5
    delta = 0.3 #  0.15-0.77
    rho = 0.2 #  0.1-0.35
    theta = 0.25 # 0.15-0.4


    nb_total_t0 = 10000
    nb_exp_t0 = 100
    nb_inf_t0 = 50
    nb_car_t0 = 50
    nb_hosp_t0 = 20
    nb_icu_t0 = 10
    nb_rec_t0 = 10
    nb_dead_t0 = 0

    # set the params required for the simulation
    params = SeirParam(tinc, tinfmild, tserint, thosp2home, thome2hosp, thosp2icu, ticu2home, tinfasy, ticu2death,
                       cont_freq, alpha, beta, delta, rho, theta,
                       nb_total_t0, nb_exp_t0, nb_car_t0, nb_inf_t0, nb_hosp_t0, nb_icu_t0, nb_rec_t0, nb_dead_t0)

    print_seir_params(params)

    # run the simulation
    result = simulate(t0=0., tmax=100., dt=0.1, params=params)

    plt.plot(result.t, result.nb_sus, label='#Suscepted')
    plt.plot(result.t, result.nb_exp, label='#Exposed')
    plt.plot(result.t, result.nb_car, label='#Carrying')
    plt.plot(result.t, result.nb_inf, label='#Infected')
    plt.plot(result.t, result.nb_hosp, label='#Hospitalzed')
    plt.plot(result.t, result.nb_icu, label='#In Intensive Care Units')
    plt.plot(result.t, result.nb_rec, label='#Recovered')
    plt.plot(result.t, result.nb_dead, label='#Died')
    plt.title("SECIR model simulation")
    plt.xlabel("Time [days]")
    plt.legend()

    plt.show()


if __name__ == "__main__":
    plot_secir()
