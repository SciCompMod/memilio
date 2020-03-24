import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


# The SIR model differential equations.
def deriv(y, t, N, alpha, beta, gamma):
    S, E, I, R = y
    # change in people susceptible to the disease
    # moderated by the number of infected people and their contact with the infected.
    dSdt = -beta * S * I / N
    # people who have been exposed to the disease
    # grows based on the contact rate and decreases based on the incubation period
    # whereby people then become infected
    dEdt = beta * S * I / N - alpha * E
    # change in infected people based on the exposed population and the incubation period
    # decreases based on the infectious period: the higher gamma is, the more quickly people die/recove
    dIdt = alpha * E - gamma * I
    # no longer infected: immune or diseased
    dRdt = gamma * I
    return dSdt, dEdt, dIdt, dRdt

# The SIR model differential equations with social distancing.
def deriv_social(y, t, N, alpha, beta, gamma, rho):
    S, E, I, R = y
    # change in people susceptible to the disease
    # moderated by the number of infected people and their contact with the infected.
    dSdt = -rho*beta * S * I / N
    # people who have been exposed to the disease
    # grows based on the contact rate and decreases based on the incubation period
    # whereby people then become infected
    dEdt = rho*beta * S * I / N - alpha * E
    # change in infected people based on the exposed population and the incubation period
    # decreases based on the infectious period: the higher gamma is, the more quickly people die/recove
    dIdt = alpha * E - gamma * I
    # no longer infected: immune or diseased
    dRdt = gamma * I
    return dSdt, dEdt, dIdt, dRdt

def plot_data(N, S,E,I,R,t):
    # Plot the data on three separate curves for S(t), I(t) and R(t)
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, fc='#dddddd', axisbelow=True)
#    ax.plot(t, S/N, 'm', alpha=0.5, lw=2, label='Susceptible')
    ax.plot(t, E, 'b', alpha=0.5, lw=2, label='Exposed')
    ax.plot(t, I, 'r', alpha=0.5, lw=2, label='Infected')
 #   ax.plot(t, R/N, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
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
    # Initial number of exposed, infected and recovered individuals, E0, I0 and R0.
    E0, I0, R0 = 1, 0, 3.5
    # Everyone else, S0, is susceptible to infection initially.
    S0 = N - I0 - R0 - E0
    # Contact rate, beta, and mean recovery rate, gamma, (in 1/days).
    #beta, gamma = 0.34, 1./10
    # End of simulation
    t_max = 200
    # time step
    dt = 0.1
    # A grid of time points (in days)
    t = np.linspace(0, t_max, int(t_max / dt))


    # A recent study of COVID-19 estimates some of these values for us (Hellewell et al. 2020),
    # so we can use some of their parameter estimates to get our model off the ground.
    # Incubation period = 5 days -> alpha = 0.2
    # This value defines how quickly the disease spreads R0 = beta/gamma
    # R0 = 3.5
    # to get 1/gamma value of 2 days, so gamma = 0.5.
    # Plugging the R0 and gamma values into Equation (6), we get an estimate of beta = 1.75.
    # https://towardsdatascience.com/social-distancing-to-slow-the-coronavirus-768292f04296
    # constraint: fixed population
    # alpha is the inverse of the incubation period (1/t_incubation)
    alpha = 0.2
    # gamma is the inverse of the mean infectious period (1/t_infectious)
    gamma = 0.5
    # beta is the average contact rate in the population
    beta = R0 * gamma
    # Coronavirus with Social Distancing
    # the term this is going to impact is our contact rate, beta.
    # a new value 0 <= rho <= 1 will capture this effect
    # 0 indicates everyone is locked down and quarantined
    # 1 is equivalent to our base case
    # hence we multiply beta with rho in our SEIR model
    rho = 1.0

    # Initial conditions vector
    y0 = S0, E0, I0, R0

    # no party
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(deriv_social, y0, t, args=(N, alpha, beta, gamma,rho))
    S, E, I, R = ret.T

    plot_data(N, S, E, I, R, t)


    # if rho is 1, we have the standard model
    # if rho == 1:
    #     # Initial conditions vector
    #     y0 = S0, E0, I0, R0
    #     # Integrate the SIR equations over the time grid, t.
    #     ret = odeint(deriv, y0, t, args=(N, alpha, beta, gamma))
    #     S, E, I, R = ret.T
    # # with social distancing effect
    # else:
    #     # Initial conditions vector
    #     y0 = S0, E0, I0, R0
    #     # Integrate the SIR equations over the time grid, t.
    #     ret = odeint(deriv_social, y0, t, args=(N, alpha, beta, gamma,rho))
    #     S, E, I, R = ret.T
    #
    # plot_data(N, S, E, I, R, t)

    # corona party
    # on day
    day_party = 30
    # with number of people
    N_party = 50
    rho_party = 0
    # the party lasts one day
    party_length = 1
    t_party = np.linspace(day_party, day_party + party_length, int(party_length / dt))

    # first everyone behaves normally
    t = np.linspace(0, day_party, int(day_party / dt))
    # Initial conditions vector
    y0 = S0, E0, I0, R0
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(deriv, y0, t, args=(N, alpha, beta, gamma))
    S, E, I, R = ret.T
    # then some people have a corona party
    # we suppose equal parts of the population join the party
    party_factor = N_party/N
    IP0 = I[-1]*party_factor
    EP0 = E[-1] * party_factor
    RP0 = I[-1] * party_factor
    SP0 = N_party - IP0 - RP0 - EP0
    print(SP0, IP0, EP0, RP0, N_party)

    # No social distancing at the party
    rho_party = 0.1
    # Initial conditions vector
    yP0 = SP0, EP0, IP0, RP0
    # Integrate the SIR equations over the time grid, t_party.
    ret = odeint(deriv_social, yP0, t_party, args=(N_party, alpha, beta, gamma,rho_party))
    SP, EP, IP, RP = ret.T
    plot_data(N_party, SP, EP, IP, RP, t_party)

    # Rest of the world behaves as always
       # Initial conditions vector
    E0 = E[-1] - EP0
    I0 = I[-1] - IP0
    R0 = R[-1] - RP0
    N_rest = N - N_party
    S0 = N_rest - I0 - R0 - E0
    print(S0, I0, E0, R0, N_rest)
    y0 = S0, E0, I0, R0
    # Integrate the SIR equations over the time grid, t_party.
    ret = odeint(deriv_social, y0, t_party, args=(N_rest, alpha, beta, gamma, rho))
    SR, ER, IR, RR = ret.T
    plot_data(N_rest, SR, ER, IR, RR, t_party)


    # now everything goes on as always
    E0 = ER[-1] + EP[-1];
    I0 = IR[-1] + IP[-1];
    R0 = RR[-1] + RP[-1];
    S0 = N - I0 - R0 - E0
    after_party = day_party+party_length
    t = np.linspace(after_party, t_max, int((t_max-after_party) / dt))
    # Initial conditions vector
    y0 = S0, E0, I0, R0
    # Integrate the SIR equations over the time grid, t.
    ret = odeint(deriv_social, y0, t, args=(N, alpha, beta, gamma, rho))
    SL, EL, IL, RL = ret.T

    t = np.linspace(0, t_max, int(t_max / dt))
    print(len(S), len(SP + SR), len(SL), len(t))
    S_final = np.concatenate((S, SP+SR, SL), axis=0) #[S, SP + SR, SL]
    E_final = np.concatenate((E, EP+ER, EL), axis=0)
    I_final = np.concatenate((I, IP+IR, IL), axis=0)
    R_final = np.concatenate((R, RP+RR, RL), axis=0)

  #S0 = N_rest - I0 - R0 - E0
    plot_data(N, S_final, E_final, I_final, R_final, t)






    plt.show()


if __name__ == "__main__":
    main()
