import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import pandas as pd

def derivative(t, X):
    '''
    A simple ODE-SEIR model

    @param X Array of values for S, E, I, and R.
    @param t Time point.
    @return Returns derivative of S, E, I, and R at t.
    '''
    S, E, I, R = X
    total_population = S + E + I + R
    derivS = - contacts * transmission_prob * S * I / total_population
    derivE = contacts * transmission_prob * S * I / total_population - E / t_e
    derivI = E / t_e - I / t_i
    derivR = I / t_i
    return np.array([derivS, derivE, derivI, derivR])

#parameter
populations = [83000]
t_e = 5.2
t_i = 6.
transmission_prob = .04
contacts = 10.

# population
E0 = 100
I0 = 50
R0 = 10
S0 = populations[0] - E0 - I0 - R0


# simulation
tmax = 100
sol = integrate.solve_ivp(derivative, (0, tmax), [S0, E0, I0, R0], atol=1e-8, rtol=1e-8, method='LSODA')

plt.figure(figsize=(10, 6))
plt.plot(sol.t, sol.y[0], label='Susceptible')
plt.plot(sol.t, sol.y[1], label='Exposed')
plt.plot(sol.t, sol.y[2], label='Infected')
plt.plot(sol.t, sol.y[3], label='Recovered')
plt.grid()
plt.legend()
plt.xlabel('Time')
plt.ylabel('Population')
plt.title('SEIR Model')
plt.show()