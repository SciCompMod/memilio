transpile = True

from integrators import euler

if transpile:
    from org.transcrypt.stubs.browser import *
    from org.transcrypt.stubs.browser import __main__, __envir__, __pragma__

# Imports for Transcrypt, skipped runtime by CPython
    if __envir__.executor_name == __envir__.transpiler_name:
        import numscrypt as np

else:

    # Imports for CPython, skipped compile time by Transcrypt
    #__pragma__ ('skip')
    import numpy as np
    from scipy.integrate import odeint
    #__pragma__ ('noskip')


# contact rate beta
# basic reproductive number of virus times inverse of mean infectious rate r_0 = 2.68 * gamma = 0.5
beta = 1.75


def deriv(y, t, N, alpha, beta, gamma):
    """
    The SIR model differential equations.
    """

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


def deriv_social(y, t, N, alpha, beta, gamma, rho):
    """
    The SEIR model differential equations with social distancing factor rho.
    with rho = 1 standard seir model (no distancing)
    can we also increase this to be > 1?
    with rho = 0 quarantine
    """
    S, E, I, R = y
    # change in people susceptible to the disease
    # moderated by the number of infectious people and their contact with the infectious.
    dSdt = -rho * beta * S * I / N
    # people who have been exposed to the disease
    # grows based on the contact rate and decreases based on the incubation period
    # whereby people then become infectious
    dEdt = rho * beta * S * I / N - alpha * E
    # change in infectious people based on the exposed population and the incubation period
    # decreases based on the infectious period: the higher gamma is, the more quickly people die/recover
    dIdt = alpha * E - gamma * I
    # no longer infected: immune or diseased
    dRdt = gamma * I
    return dSdt, dEdt, dIdt, dRdt


def my_odeint(func, y0, t, args=()):
    """
    Alternative to scipy odeint
    """

    t_min = t[0]
    t_max = t[-1]
    dt = t[1] - t[0]

    def fun_helper(t, y):
        return func(y, t, *args)

    tp, vals = euler(fun_helper, y0, (t_min, t_max), dt)

    return vals[:, 0:t.shape[0]].T

def seir_solve(N, S0, E0, I0, R0, t_max, dt, alpha, gamma, rho, t_min=0):
    # Initial conditions vector
    y0 = S0, E0, I0, R0
    # A grid of time points (in days)

    t = np.linspace(t_min, t_max, int((t_max - t_min) / dt))
    # Integrate the SEIR equations over the time grid, t.
    # ret2 = odeint(deriv_social, y0, t, args=(N, alpha, beta, gamma, rho))
    ret = my_odeint(deriv_social, y0, t, args=(N, alpha, beta, gamma, rho))

    S, E, I, R = ret.T
    # plot_data(N, S, E, I, R, t)
    return S, E, I, R


class EpiParms:
    # Define parameters
    # Total population, N.
    N = 10000
    # Initial number of exposed, infectious and recovered individuals, E0, I0 and R0.
    E0, I0, R0 = 1, 0, 0

    # End of simulation
    t_max = 200
    # time step
    dt = 0.1
    # Assume an incubation period of 5.2 days
    t_incubation = 5.2
    # Assume infectious period of 2 days
    t_infectious = 2

    # Parameter for Social Distancing
    # a new value 0 <= rho <= 1 will capture this effect
    # the term this is going to impact is our contact rate, beta.
    # 0 indicates everyone is locked down and quarantined
    # 1 is equivalent to our base case
    # hence we multiply beta with rho in our SEIR model
    rho = 0.5


def epi(parameters: EpiParms):
    # Define parameters
    N = parameters.N
    E0 = parameters.E0
    I0 = parameters.I0
    R0 = parameters.R0
    t_max = parameters.t_max
    dt = parameters.dt
    rho = parameters.rho

    # Everyone else, S0, is susceptible to infection initially.
    # Constraint: fixed population
    S0 = parameters.N - parameters.I0 - parameters.R0 - parameters.E0

    # alpha is the inverse of the incubation period (1/t_incubation)
    alpha = 1.0 / parameters.t_incubation
    # gamma is the mean recovery rate
    # inverse of the mean infectious period (1/t_infectious)
    gamma = 1.0 / parameters.t_infectious

    # solve the SEIR model and plot the results
    SA, EA, IA, RA = seir_solve(N, S0, E0, I0, R0, t_max, dt, alpha, gamma, rho, t_min=0)

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
    S_till_party, E_till_party, I_till_party, R_till_party = seir_solve(N, S0, E0, I0, R0, day_party, dt, alpha, gamma,
                                                                        rho, t_min=0)
    # print("before party", S_till_party[-1], I_till_party[-1], E_till_party[-1], R_till_party[-1], N)

    # people have a corona party
    # we suppose equal parts of the population join the party
    party_factor = N_party / N
    # we get these equal parts from the solution up until this day
    IP0 = I_till_party[-1] * party_factor
    EP0 = E_till_party[-1] * party_factor
    RP0 = R_till_party[-1] * party_factor
    SP0 = N_party - IP0 - RP0 - EP0
    # corona party solve
    S_party, E_party, I_party, R_party = seir_solve(N_party, SP0, EP0, IP0, RP0, day_party + length_party, dt, alpha,
                                                    gamma, rho_party, t_min=day_party)
    # print("after party", S_party[-1], I_party[-1], E_party[-1], R_party[-1], N_party)

    # Rest of the world behaves as always with more social distance
    # We deduce the people attending the party
    rest_factor = (N - N_party) / N
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
    # print("all after party", S0, "=", SR[-1]+S_party[-1], E0,I0,S0, N, ER[-1], E_party[-1])
    SL, EL, IL, RL = seir_solve(N, S0, E0, I0, R0, t_max, dt, alpha, gamma, rho,
                                t_min=day_party + length_party)

    # now put everything together
    t = np.linspace(0, t_max, int(t_max / dt))
    S_final = np.concatenate((S_till_party, S_party + SR, SL), axis=0)
    E_final = np.concatenate((E_till_party, E_party + ER, EL), axis=0)
    I_final = np.concatenate((I_till_party, I_party + IR, IL), axis=0)
    R_final = np.concatenate((R_till_party, R_party + RR, RL), axis=0)

    E_party_max = np.amax(E_final)
    I_party_max = np.amax(I_final)
    R_party_max = np.amax(R_final)
    EA_max = np.amax(EA)
    IA_max = np.amax(IA)
    RA_max = np.amax(RA)
    # people who all in all had contact with Corona
    r1 = R_final[-1] - RA[-1]
    print("Population:", N, "Day of party:", day_party, "Party people:", N_party, "Death rate assumption: 2%")
    print("At the worst day in our scenario, ", int(E_party_max + I_party_max), " people are infected.")
    print("A the end of the crisis", int(R_party_max), "people have been infected,")
    print("without the party it would have been", round(r1, 3), "less.")
    print("You think that this does not make a difference?")
    # Population of Germany 83 783 942
    factor = 83783942 / N
    print("In all of Germany this makes a difference of", int(r1 * factor), "people who are infected.")
    print(int(r1 * factor * 2 / 100), "people will have died due to the parties.")

    return S_final, E_final, I_final, R_final, t
