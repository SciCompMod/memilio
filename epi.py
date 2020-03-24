import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt


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
t_max = 100
# time step
dt = 0.1
# A grid of time points (in days)
t = np.linspace(0, t_max, int(t_max / dt) + 1)


# one exposed person S_0 = 1–1/N, E_0 = 1/N, I_0 = 0, R_0 = 0
init_vals = 1 - 1 / N, 1 / N, 0, 0

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



#params = alpha, beta, gamma
# The SIR model differential equations.
def deriv(y, t, N, alpha, beta, gamma):
    S, E, I, R = y
    # change in people susceptible to the disease
    # moderated by the number of infected people and their contact with the infected.
    dSdt = -beta * S * I / N
    # people who have been exposed to the disease
    # grows based on the contact rate and decreases based on the incubation period
    # whereby people then become infected
    dEdt = beta * S*I/N - alpha*E
    # change in infected people based on the exposed population and the incubation period
    # decreases based on the infectious period: the higher gamma is, the more quickly people die/recove
    dIdt = alpha*E - gamma * I
    # no longer infected: immune or diseased
    dRdt = gamma * I
    return dSdt, dEdt, dIdt, dRdt

# Initial conditions vector
y0 = S0, E0, I0, R0
# Integrate the SIR equations over the time grid, t.
ret = odeint(deriv, y0, t, args=(N, alpha, beta, gamma))
S, E, I, R = ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, fc='#dddddd', axisbelow=True)
#ax.plot(t, S/N, 'm', alpha=0.5, lw=2, label='Susceptible')
ax.plot(t, E/N, 'b', alpha=0.5, lw=2, label='Exposed')
ax.plot(t, I/N, 'r', alpha=0.5, lw=2, label='Infected')
#ax.plot(t, R/N, 'g', alpha=0.5, lw=2, label='Recovered with immunity')
ax.set_xlabel('Days')
ax.set_ylabel('Population Fraction ')
#ax.set_ylim(bottom,top)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('top', 'right', 'bottom', 'left'):
    ax.spines[spine].set_visible(False)
plt.show()

