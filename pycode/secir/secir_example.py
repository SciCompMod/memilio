from epidemiology.secir import Damping, SecirParams, print_secir_params, simulate, StageTimes, Probabilities, Populations
from matplotlib import pyplot as plt


def plot_secir():
    """
    Simulates the secir model using the c++ code
    and plots the results
    """

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
    times.set_cont_freq(0.5)  # 0.2-0.75

    probs = Probabilities()
    probs.set_asymp_per_infectious(0.09)  # 0.01-0.16
    probs.set_risk_from_symptomatic(0.25)  # 0.05-0.5
    probs.set_hospitalized_per_infectious(0.2)  # 0.1-0.35
    probs.set_icu_per_hospitalized(0.25)  # 0.15-0.4
    probs.set_dead_per_icu(0.3)  # 0.15-0.77

    people = Populations()
    people.set_total_t0(10000)
    people.set_exposed_t0(100)
    people.set_carrier_t0(50)
    people.set_infectious_t0(50)
    people.set_hospital_t0(20)
    people.set_icu_t0(10)
    people.set_recovered_t0(10)
    people.set_dead_t0(0)

    # set the params required for the simulation
    params = SecirParams()
    params.times = times
    params.probabilities = probs
    params.populations = people

    # emulate some mitigations
    params.dampings.add(Damping(23., 0.8))
    params.dampings.add(Damping(25., 0.75))
    params.dampings.add(Damping(27., 0.7))

    print_secir_params(params)

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
