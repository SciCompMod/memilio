from epi import epi, EpiParms
import matplotlib.pyplot as plt


def plot_data(N, S, E, I, R, t):
    """
    Plot the data on three separate curves for S(t), I(t) and R(t)
    """
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
    parameters = EpiParms()

    # Define parameters
    # Total population, N.
    parameters.N = 10000
    # Initial number of exposed, infectious and recovered individuals, E0, I0 and R0.
    parameters.E0, parameters.I0, parameters.R0 = 1, 0, 0

    # End of simulation
    parameters.t_max = 200
    # time step
    parameters.dt = 0.1
    # Assume an incubation period of 5.2 days
    parameters.t_incubation = 5.2
    # Assume infectious period of 2 days
    parameters.t_infectious = 2

    # Parameter for Social Distancing
    # a new value 0 <= rho <= 1 will capture this effect
    # the term this is going to impact is our contact rate, beta.
    # 0 indicates everyone is locked down and quarantined
    # 1 is equivalent to our base case
    # hence we multiply beta with rho in our SEIR model
    parameters.rho = 0.5

    # perform the simulation
    S_final, E_final, I_final, R_final, t = epi(parameters)

    plot_data(parameters.N, S_final, E_final, I_final, R_final, t)
    plt.show()


if __name__ == "__main__":
    main()
