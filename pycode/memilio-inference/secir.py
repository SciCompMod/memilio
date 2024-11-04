import numpy as np
import numpy.typing as npt
from scipy import stats
from enum import Enum
from functools import partial

from memilio.simulation import AgeGroup, Damping, RungeKuttaCashKarp54IntegratorCore, RKIntegratorCore
from memilio.simulation.osir import (InfectionState, Model, interpolate_simulation_result,
                                     simulate_flows)
from prior import (UnboundParameter, LambdaParameter, DelayParameter,
                   InterventionChangePointParameter, WeeklyModulationAmplitudeParameter,
                   ScaleMultiplicativeReportingNoiseParameter, ModelStrategy)

import matplotlib.pyplot as plt

alpha_f = (0.7**2)*((1-0.7)/(0.17**2) - (1-0.7))
beta_f = alpha_f*(1/0.7 - 1)


class SECIRParameterNames(Enum):
    LAMBDA_0 = r'$\lambda_0$'
    MU = r'$\mu$'
    E0 = r'$E_0$'
    ALPHA = r'$\alpha$'
    BETA = r'$\beta$'
    GAMMA = r'$\gamma$'
    ETA = r'$\eta$'
    THETA = r'$\theta$'
    DELTA = r'$\delta$'
    D = r'$\d$'
    T1 = r'$t_1$'
    T2 = r'$t_2$'
    T3 = r'$t_3$'
    T4 = r'$t_4$'
    LAMBDA_1 = r'$\lambda_1$'
    LAMBDA_2 = r'$\lambda_2$'
    LAMBDA_3 = r'$\lambda_3$'
    LAMBDA_4 = r'$\lambda_4$'
    F_I = r'$f_i$'
    PHI_I = r'$\phi_i$'
    DELAY_I = r'$D_i$'
    PSI_I = r'$\psi_i$'
    F_D = r'$f_d$'
    PHI_D = r'$\phi_i$'
    DELAY_D = r'$D_d$'
    PSI_D = r'$\psi_d$'


class SECIRStrategy(ModelStrategy):
    @staticmethod
    # Possible option to draw without redraw
    def add_base(prior_array: list[UnboundParameter]) -> None:
        # Params same as SIR
        prior_array.append(LambdaParameter(distribution=partial(
            np.random.lognormal, mean=np.log(1.2), sigma=0.5), name=SECIRParameterNames.LAMBDA_0.value))
        prior_array.append(UnboundParameter(distribution=partial(
            np.random.lognormal, mean=np.log(8), sigma=0.2), name=SECIRParameterNames.MU.value))
        prior_array.append(UnboundParameter(distribution=partial(
            np.random.gamma, shape=2, scale=30), name=SECIRParameterNames.E0.value))

        # New Params
        prior_array.append(UnboundParameter(distribution=partial(
            np.random.uniform, low=0.005, high=0.99), name=SECIRParameterNames.ALPHA.value))
        prior_array.append(UnboundParameter(distribution=partial(
            np.random.lognormal, mean=np.log(0.25), sigma=0.3), name=SECIRParameterNames.BETA.value))
        prior_array.append(UnboundParameter(distribution=partial(
            np.random.lognormal, mean=np.log(1/6.5), sigma=0.5), name=SECIRParameterNames.GAMMA.value))
        prior_array.append(UnboundParameter(distribution=partial(
            np.random.lognormal, mean=np.log(1/3.2), sigma=0.3), name=SECIRParameterNames.ETA.value))
        prior_array.append(UnboundParameter(distribution=partial(
            np.random.uniform, low=1/14, high=1/3), name=SECIRParameterNames.THETA.value))
        prior_array.append(UnboundParameter(distribution=partial(
            np.random.uniform, low=0.01, high=0.3), name=SECIRParameterNames.DELTA.value))
        prior_array.append(UnboundParameter(distribution=partial(
            np.random.uniform, low=1/14, high=1/3), name=SECIRParameterNames.D.value))

    @staticmethod
    # Possible option to draw without redraw
    def add_intervention(prior_array: list[UnboundParameter]) -> None:
        # Params same as SIR
        prior_array.append(InterventionChangePointParameter(distribution=partial(
            np.random.normal, loc=8, scale=3), name=SECIRParameterNames.T1.value))
        prior_array.append(InterventionChangePointParameter(distribution=partial(
            np.random.normal, loc=15, scale=1), name=SECIRParameterNames.T2.value))
        prior_array.append(InterventionChangePointParameter(distribution=partial(
            np.random.normal, loc=22, scale=1), name=SECIRParameterNames.T3.value))
        prior_array.append(InterventionChangePointParameter(distribution=partial(
            np.random.normal, loc=66, scale=1), name=SECIRParameterNames.T4.value))
        prior_array.append(LambdaParameter(distribution=partial(
            np.random.lognormal, mean=np.log(0.6), sigma=0.5), name=SECIRParameterNames.LAMBDA_1.value))
        prior_array.append(LambdaParameter(distribution=partial(
            np.random.lognormal, mean=np.log(0.3), sigma=0.5), name=SECIRParameterNames.LAMBDA_2.value))
        prior_array.append(LambdaParameter(distribution=partial(
            np.random.lognormal, mean=np.log(0.1), sigma=0.5), name=SECIRParameterNames.LAMBDA_3.value))
        prior_array.append(LambdaParameter(distribution=partial(
            np.random.lognormal, mean=np.log(0.1), sigma=0.5), name=SECIRParameterNames.LAMBDA_4.value))

    @staticmethod
    # Possible option to draw without redraw
    def add_observation(prior_array: list[UnboundParameter], sim_diff: int) -> None:
        prior_array.append(WeeklyModulationAmplitudeParameter(distribution=partial(
            np.random.beta, a=alpha_f, b=beta_f), name=SECIRParameterNames.F_I.value))
        prior_array.append(UnboundParameter(
            distribution=stats.vonmises(kappa=0.01).rvs, name=SECIRParameterNames.PHI_I.value))
        prior_array.append(DelayParameter(distribution=partial(
            np.random.lognormal, mean=np.log(8), sigma=0.2), name=SECIRParameterNames.DELAY_I.value, sim_diff=sim_diff))
        prior_array.append(ScaleMultiplicativeReportingNoiseParameter(distribution=partial(
            np.random.gamma, shape=1, scale=5), name=SECIRParameterNames.PSI_I.value))

        prior_array.append(WeeklyModulationAmplitudeParameter(distribution=partial(
            np.random.beta, a=alpha_f, b=beta_f), name=SECIRParameterNames.F_D.value))
        prior_array.append(UnboundParameter(
            distribution=stats.vonmises(kappa=0.01).rvs, name=SECIRParameterNames.PHI_D.value))
        prior_array.append(DelayParameter(distribution=partial(
            np.random.lognormal, mean=np.log(8), sigma=0.2), name=SECIRParameterNames.DELAY_D.value, sim_diff=sim_diff))
        prior_array.append(ScaleMultiplicativeReportingNoiseParameter(distribution=partial(
            np.random.gamma, shape=1, scale=5), name=SECIRParameterNames.PSI_D.value))


def stationary_SECIR(params: list[float], N: int, T: int, intervention_model: bool, observation_model: bool, param_names: list[str], fixed_params: dict[str, float] = dict(), sim_diff: int = 16) -> npt.NDArray[np.float64]:
    """Performs a forward simulation from the stationary SECIR model given a random draw from the prior."""

    # Convert params array to a dictionary using param_names
    dynamic_params = {name: params[i] for i, name in enumerate(param_names)}

    # Merge with fixed parameters, overriding dynamic params where fixed
    combined_params = {**dynamic_params, **fixed_params}

    # Extract parameters using the Enum
    lambd0 = combined_params[SECIRParameterNames.LAMBDA_0.value]
    mu = combined_params[SECIRParameterNames.MU.value]
    alpha = combined_params[SECIRParameterNames.ALPHA.value]
    beta = combined_params[SECIRParameterNames.BETA.value]
    gamma = combined_params[SECIRParameterNames.GAMMA.value]
    eta = combined_params[SECIRParameterNames.ETA.value]
    theta = combined_params[SECIRParameterNames.THETA.value]
    delta = combined_params[SECIRParameterNames.DELTA.value]
    d = combined_params[SECIRParameterNames.D.value]
    E0 = combined_params[SECIRParameterNames.E0.value]
    # Should IO also be a constraint instead of set to bouandry?
    E0 = max(1, np.round(E0))

    if intervention_model:
        t1 = combined_params[SECIRParameterNames.T1.value]
        t2 = combined_params[SECIRParameterNames.T2.value]
        t3 = combined_params[SECIRParameterNames.T3.value]
        t4 = combined_params[SECIRParameterNames.T4.value]
        # Round integer parameters
        t1, t2, t3, t4 = int(round(t1)), int(
            round(t2)), int(round(t3)), int(round(t4))

        lambd1 = combined_params[SECIRParameterNames.LAMBDA_1.value]
        lambd2 = combined_params[SECIRParameterNames.LAMBDA_2.value]
        lambd3 = combined_params[SECIRParameterNames.LAMBDA_3.value]
        lambd4 = combined_params[SECIRParameterNames.LAMBDA_4.value]

    if observation_model:
        f_i = combined_params[SECIRParameterNames.F_I.value]
        phi_i = combined_params[SECIRParameterNames.PHI_I.value]
        D_i = combined_params[SECIRParameterNames.DELAY_I.value]
        scale_I = combined_params[SECIRParameterNames.PSI_I.value]
        D_i = int(round(D_i))

        f_d = combined_params[SECIRParameterNames.F_D.value]
        phi_d = combined_params[SECIRParameterNames.PHI_D.value]
        D_d = combined_params[SECIRParameterNames.DELAY_D.value]
        scale_d = combined_params[SECIRParameterNames.PSI_D.value]
        D_d = int(round(D_d))

    # Maybe check for missing parameter

    # Configure model
    num_groups = 1
    secir_model = Model(num_groups)
    A0 = AgeGroup(0)

    # Calculate lambda arrays
    # Lambda0 is the initial contact rate which will be consecutively
    # reduced via the government measures
    sim_lag = sim_diff - 1

    # Initial conditions
    secir_model.populations[A0, InfectionState.Exposed] = E0
    secir_model.populations[A0, InfectionState.Infected] = 0
    secir_model.populations[A0, InfectionState.Recovered] = 0
    secir_model.populations[A0, InfectionState.Dead] = 0
    secir_model.populations.set_difference_from_total(
        (A0, InfectionState.Susceptible), N)

    # Initialize Parameters
    secir_model.parameters.TimeInfected[A0] = mu
    secir_model.parameters.TransmissionProbabilityOnContact[A0] = 1
    # Very weird way of setting dampings and differnt from paper, because of no time till NPIs fully take effect
    # Also damping could increase the beta value, should that be possible?

    def calc_damping_value_from_lambda(_lambdt: float, _lambd0: float = 1, _baseline: float = 1, _minimum: float = 0) -> float:
        # cf = bl - (dampingvalue * (bl - min))
        # -> bl - cf = (dampingvalue * (bl - min)
        # -> (bl - cf) / (bl - min) = dampingvalue
        # lambdt = cf * lambd0 -> cf = lambdt / lambd0
        # => dampingvalue = (bl - (lambdt / lambd0)) / (bl - min)
        return (_baseline - (_lambdt / _lambd0)) / (_baseline - _minimum)

    secir_model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
        (num_groups, num_groups))
    secir_model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.zeros(
        (num_groups, num_groups))
    secir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
        coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd0), t=0, level=0, type=0))

    if intervention_model:
        # delta_t1, delta_t2, delta_t3, delta_t4 = int(round(delta_t1)), int(round(delta_t2)), int(round(delta_t3)), int(round(delta_t4))
        secir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
            coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd1), t=t1+sim_lag, level=0, type=0))
        secir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
            coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd2), t=t2+sim_lag, level=0, type=0))
        secir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
            coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd3), t=t3+sim_lag, level=0, type=0))
        secir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
            coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd4), t=t4+sim_lag, level=0, type=0))

    # Check logical constraints to parameters, should we really use apply_constraints?!
    secir_model.apply_constraints()

    # Reported new cases
    I_data = np.zeros(T)
    fs_i = np.zeros(T)

    # Run Simulation
    integrator = RKIntegratorCore(dt_max=1)

    (result, flows) = simulate_flows(
        0, T + sim_lag - 1, 1, secir_model, integrator)

    # (result, flows) = simulate_flows(
    #     0, T + sim_lag - 1, 1, secir_model)
    # interpolate results
    flows = interpolate_simulation_result(flows)

    if observation_model:

        timepoints = flows.get_num_time_points()

        # Adding new cases with delay D
        # Note, we assume the same delay
        I_data = np.diff(flows.as_ndarray()[
                         1, (sim_lag-D_i-1):(sim_lag-D_i)+T])
        I_data = np.clip(I_data, 10 ** -14, N)

        # Compute lags
        fs_i = (1-f_i)*(1 -
                        np.abs(np.sin((np.pi/7) * np.arange(0, timepoints-sim_lag, 1) - 0.5*phi_i)))

        # Compute weekly modulation
        I_data = (1-fs_i) * I_data

        # check for negative values
        try:
            scale = np.sqrt(I_data)*scale_I
            assert np.all(scale >= 0)
        except AssertionError as e:
            print('Invalid value simulated...return nan')
            return np.stack(([np.nan] * 81, )).T

        # Add noise
        I_data = stats.t(df=4, loc=I_data, scale=np.sqrt(I_data)*scale_I).rvs()
    else:
        I_data = np.diff(flows.as_ndarray()[1, :T+1])
    # bound all negative values to 0
    I_data = np.clip(I_data, 10 ** -14, N)
    return np.stack((I_data, )).T
    # at this point code example uses cum sum instead of new
    # return np.stack((flows[:,0], )).T # np.stack((result[:,1])).T also wrong because its not cumsum


if __name__ == "__main__":

    results = simulator_SIR([0.7, 4, 200], 83e6, 81, False, False, [
        SECIRParameterNames.LAMBDA_0.value, SECIRParameterNames.MU.value, SECIRParameterNames.I0.value])

    plt.plot(results)
    plt.savefig("simulator_SIR.png")
