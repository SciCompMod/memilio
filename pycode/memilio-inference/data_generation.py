import numpy as np
from memilio.simulation import AgeGroup, Damping
from memilio.simulation.osir import (InfectionState, Model, interpolate_simulation_result,
                                     simulate_flows)


def stationary_SIR(params, N, T, observation_model, sim_diff=16):
    """Performs a forward simulation from the stationary SIR model given a random draw from the prior."""

    # Extract parameters
    lambd0, mu, I0, t1, t2, t3, t4, lambd1, lambd2, lambd3, lambd4, f_i, phi_i, D_i, scale_I = params

    # Round integer parameters
    t1, t2, t3, t4 = int(round(t1)), int(
        round(t2)), int(round(t3)), int(round(t4))
    # delta_t1, delta_t2, delta_t3, delta_t4 = int(round(delta_t1)), int(round(delta_t2)), int(round(delta_t3)), int(round(delta_t4))
    # Should IO also be a constraint instead of set to bouandry?
    I0 = max(1, np.round(I0))
    D_i = int(round(D_i))

    # Configure model
    num_groups = 1
    sir_model = Model(num_groups)
    A0 = AgeGroup(0)

    # Calculate lambda arrays
    # Lambda0 is the initial contact rate which will be consecutively
    # reduced via the government measures
    sim_lag = sim_diff - 1
    # lambd_arr = calc_lambda_array(sim_lag, lambd0, lambd1, lambd2, lambd3, lambd4,
    #                               t1, t2, t3, t4, delta_t1, delta_t2, delta_t3, delta_t4, T)

    # Initial conditions
    sir_model.populations[A0, InfectionState.Infected] = I0
    sir_model.populations[A0, InfectionState.Recovered] = 0
    sir_model.populations.set_difference_from_total(
        (A0, InfectionState.Susceptible), N)

    # Initialize Parameters
    sir_model.parameters.TimeInfected[A0] = mu
    sir_model.parameters.TransmissionProbabilityOnContact[A0] = 1
    # Very weird way of setting dampings and differnt from paper, because of no time till NPIs fully take effect
    # Also damping could increase the beta value, should that be possible?

    def calc_damping_value_from_lambda(_lambdt, _lambd0=1, _baseline=1, _minimum=0):
        # cf = bl - (dampingvalue * (bl - min))
        # -> bl - cf = (dampingvalue * (bl - min)
        # -> (bl - cf) / (bl - min) = dampingvalue
        # lambdt = cf * lambd0 -> cf = lambdt / lambd0
        # => dampingvalue = (bl - (lambdt / lambd0)) / (bl - min)
        return (_baseline - (_lambdt / _lambd0)) / (_baseline - _minimum)

    sir_model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
        (num_groups, num_groups))
    sir_model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.zeros(
        (num_groups, num_groups))
    sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd0), t=0, level=0, type=0))
    sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd1), t=t1+sim_lag, level=0, type=0))
    sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd2), t=t2+sim_lag, level=0, type=0))
    sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd3), t=t3+sim_lag, level=0, type=0))
    sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd4), t=t4+sim_lag, level=0, type=0))

    # Check logical constraints to parameters
    sir_model.apply_constraints()

    # Reported new cases
    I_data = np.zeros(T)
    fs_i = np.zeros(T)

    # Run Simulation, maybe i need a flow model?
    (result, flows) = simulate_flows(0, T + sim_lag - 1, 1, sir_model)
    # interpolate results
    flows = interpolate_simulation_result(flows)

    if observation_model:

        # Calculate delay and changepoints !vectorice this!
        for t in range(flows.get_num_time_points()):

            # From here, start adding new cases with delay D
            # Note, we assume the same delay
            if t >= sim_lag:

                # Compute lags and add to data arrays
                fs_i[t-sim_lag] = (1-f_i)*(1 -
                                           np.abs(np.sin((np.pi/7) * (t-sim_lag) - 0.5*phi_i)))

        I_data = np.diff(flows.as_ndarray()[
                         1, (sim_lag-D_i-1):(sim_lag-D_i)+T])
        assert (fs_i == fs_i).all()

        # Compute weekly modulation
        I_data = (1-fs_i) * I_data

        # check for negative values
        try:
            scale = np.sqrt(I_data)*scale_I
            assert np.all(scale > 0)
        except AssertionError as e:
            print("Invalid value simulated...return nan")
            return np.stack(([np.nan] * 81, )).T

        # Add noise
        # I_data = stats.t(df=4, loc=I_data, scale=np.sqrt(I_data)*scale_I).rvs()

        # bound all negative values to 0
        I_data = np.clip(I_data, 0, N)
    else:
        I_data = np.diff(flows.as_ndarray()[1, :T+1])
    return np.stack((I_data, )).T
    # at this point code example uses cum sum instead of new
    # return np.stack((flows[:,0], )).T # np.stack((result[:,1])).T also wrong because its not cumsum
