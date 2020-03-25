import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
# contact rate beta
# basic reproductive number of virus times inverse of mean infectious rate r_0 = 2.68 * gamma = 0.5
beta = 1.75
# The SIR model differential equations.
def deriv(y, t, N, alpha, beta, gamma):
    S, E, I, R = y
    # change in people susceptible to the disease
    # moderated by the number of infectious people and their contact with the infectious.
    dSdt = -beta * S * I / N
    # people who have been exposed to the disease
    # grows based on the contact rate and decreases based on the incubation period
    # whereby people then become infectious
    dEdt = beta * S * I / N - alpha * E
    # change in infectious people based on the exposed population and the incubation period
    # decreases based on the infectious period: the higher gamma is, the more quickly people die/recover
    dIdt = alpha * E - gamma * I
    # no longer infected: immune or diseased
    dRdt = gamma * I
    return dSdt, dEdt, dIdt, dRdt

# The SEIR model differential equations with social distancing factor rho.
# with rho = 1 standard seir model (no distancing)
# can we also increase this to be > 1?
# with rho = 0 quarantine
def deriv_social(y, t, N, alpha, beta, gamma, rho):
    S, E, I, R = y
    # change in people susceptible to the disease
    # moderated by the number of infectious people and their contact with the infectious.
    dSdt = -rho*beta * S * I / N
    # people who have been exposed to the disease
    # grows based on the contact rate and decreases based on the incubation period
    # whereby people then become infectious
    dEdt = rho*beta * S * I / N - alpha * E
    # change in infectious people based on the exposed population and the incubation period
    # decreases based on the infectious period: the higher gamma is, the more quickly people die/recover
    dIdt = alpha * E - gamma * I
    # no longer infected: immune or diseased
    dRdt = gamma * I
    return dSdt, dEdt, dIdt, dRdt

def seir_solve(N, S0, E0, I0, R0, t_max, dt, alpha, gamma, rho, t_min=0):
    # Initial conditions vector
    y0 = S0, E0, I0, R0
    # A grid of time points (in days)
    t = np.linspace(t_min, t_max, int((t_max-t_min) / dt))
    # Integrate the SEIR equations over the time grid, t.
    ret = odeint(deriv_social, y0, t, args=(N, alpha, beta, gamma, rho))
    S, E, I, R = ret.T
    # plot_data(N, S, E, I, R, t)
    return S, E, I, R

def plot_data(N, S, E, I, R, t):
    # Plot the data on three separate curves for S(t), I(t) and R(t)
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, fc='#dddddd', axisbelow=True)
    ax.plot(t, S+E+I+R, 'm', alpha=0.5, lw=2, label='S+E+I+R')
    ax.plot(t, S, 'g', alpha=0.5, lw=2, label='Susceptible')
    ax.plot(t, E, 'c', alpha=0.5, lw=2, label='Exposed')
    ax.plot(t, I, 'r', alpha=0.5, lw=2, label='Infectious')
    ax.plot(t, R-R*2/100, 'k', alpha=0.5, lw=2, label='Recovered')
    ax.plot(t, R*2/100, 'b', alpha=0.5, lw=2, label='Dead')
    ax.set_xlabel('Days')
    ax.set_ylabel('Population')
    #ax.set_ylim(bottom,top)
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('top', 'right', 'bottom', 'left'):
        ax.spines[spine].set_visible(False)
    plt.draw()

def main():
    # Define parameters
    # Total population, N.
    N = 10000
    # Initial number of exposed, infectious and recovered individuals, E0, I0 and R0.
    E0, I0, R0 = 1, 0, 0
    # Everyone else, S0, is susceptible to infection initially.
    # Constraint: fixed population
    S0 = N - I0 - R0 - E0
    # End of simulation
    t_max = 200
    # time step
    dt = 0.1
    # Assume an incubation period of 5.2 days
    t_incubation = 5.2
    # Assume infectious period of 2 days
    t_infectious = 2
    # alpha is the inverse of the incubation period (1/t_incubation)
    alpha = 1.0 / t_incubation
    # gamma is the mean recovery rate
    # inverse of the mean infectious period (1/t_infectious)
    gamma = 1.0 / t_infectious
    # Parameter for Social Distancing
    # a new value 0 <= rho <= 1 will capture this effect
    # the term this is going to impact is our contact rate, beta.
    # 0 indicates everyone is locked down and quarantined
    # 1 is equivalent to our base case
    # hence we multiply beta with rho in our SEIR model
    rho = 0.5

    # solve the SEIR model and plot the results
    SA,EA,IA,RA = seir_solve(N, S0, E0, I0, R0, t_max, dt, alpha, gamma, rho, t_min=0)

    # corona party params
    # on the day
    day_party = 100
    # with number of people
    N_party = 50
    # with little social distancing
    rho_party = 1.0
    # the party lasts all night
    length_party = 0.5

    # solve with whole population until day of party
    S_till_party, E_till_party, I_till_party, R_till_party = seir_solve(N, S0, E0, I0, R0, day_party, dt, alpha, gamma, rho, t_min=0)
    # print("before party", S_till_party[-1], I_till_party[-1], E_till_party[-1], R_till_party[-1], N)

    # people have a corona party
    # we suppose equal parts of the population join the party
    party_factor = N_party/N
    # we get these equal parts from the solution up until this day
    IP0 = I_till_party[-1] * party_factor
    EP0 = E_till_party[-1] * party_factor
    RP0 = R_till_party[-1] * party_factor
    SP0 = N_party - IP0 - RP0 - EP0
    # corona party solve
    S_party, E_party, I_party, R_party = seir_solve(N_party, SP0, EP0, IP0, RP0, day_party + length_party, dt, alpha, gamma, rho_party, t_min=day_party)
    # print("after party", S_party[-1], I_party[-1], E_party[-1], R_party[-1], N_party)


    # Rest of the world behaves as always with more social distance
    # We deduce the people attending the party
    rest_factor = (N-N_party)/N
    I0 = I_till_party[-1] * rest_factor
    E0 = E_till_party[-1] * rest_factor
    R0 = R_till_party[-1] * rest_factor
    N_rest = N - N_party
    S0 = N_rest - I0 - R0 - E0
    # print("after party rest",S_till_party[-1] * rest_factor, "=", S0, I0, E0, R0, N_rest)
    # Integrate the SEIR equations over the time grid t_party.
    SR, ER, IR, RR = seir_solve(N_rest, S0, E0, I0, R0, day_party + length_party, dt, alpha, gamma, rho,
                                t_min=day_party)

    # now everything goes on as always
    # people attending the party and rest go together once more
    E0 = ER[-1] + E_party[-1];
    I0 = IR[-1] + I_party[-1];
    R0 = RR[-1] + R_party[-1];
    S0 = N - I0 - R0 - E0
    #print("all after party", S0, "=", SR[-1]+S_party[-1], E0,I0,S0, N, ER[-1], E_party[-1])
    SL, EL, IL, RL = seir_solve(N, S0, E0, I0, R0, t_max, dt, alpha, gamma, rho,
                                t_min=day_party+length_party)


    # now put everything together
    t = np.linspace(0, t_max, int(t_max / dt))
    S_final = np.concatenate((S_till_party, S_party+SR, SL), axis=0)
    E_final = np.concatenate((E_till_party, E_party+ER, EL), axis=0)
    I_final = np.concatenate((I_till_party, I_party+IR, IL), axis=0)
    R_final = np.concatenate((R_till_party, R_party+RR, RL), axis=0)

    plot_data(N, S_final, E_final, I_final, R_final, t)
    plt.show()

    E_party_max = np.amax(E_final)
    I_party_max = np.amax(I_final)
    R_party_max = np.amax(R_final)
    EA_max = np.amax(EA)
    IA_max = np.amax(IA)
    RA_max = np.amax(RA)
    # people who all in all had contact with Corona
    r1 = R_final[-1]-RA[-1]
    print("Population:", N, "Day of party:", day_party, "Party people:", N_party, "Death rate assumption: 2%")
    print("At the worst day in our scenario, ", int(E_party_max + I_party_max), " people are infected.")
    print("A the end of the crisis", int(R_party_max), "people have been infected,")
    print("without the party it would have been",  round(r1,3), "less.")
    print("You think that this does not make a difference?")
    # Population of Germany 83 783 942
    factor = 83783942/N
    print("In all of Germany this makes a difference of", int(r1 * factor), "people who are infected.")
    print(int(r1*factor*2/100), "people will have died due to the parties.")
    exit(0)


if __name__ == "__main__":
    main()
