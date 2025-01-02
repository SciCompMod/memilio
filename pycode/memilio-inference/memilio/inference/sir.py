import numpy as np
import numpy.typing as npt
from scipy import stats
from enum import Enum
from functools import partial

from memilio.simulation import AgeGroup, Damping, RungeKuttaCashKarp54IntegratorCore, RKIntegratorCore
from memilio.simulation.osir import (InfectionState, Model, interpolate_simulation_result,
                                     simulate_flows)
from memilio.inference.prior import (UnboundParameter, LambdaParameter, DelayParameter,
                                     InterventionChangePointParameter, WeeklyModulationAmplitudeParameter,
                                     ScaleMultiplicativeReportingNoiseParameter, ModelStrategy)

import matplotlib.pyplot as plt

alpha_f = (0.7**2)*((1-0.7)/(0.17**2) - (1-0.7))
beta_f = alpha_f*(1/0.7 - 1)


# Vorberechnung für Beta-Verteilung
alpha_f = (0.7**2) * ((1 - 0.7) / (0.17**2) - (1 - 0.7))
beta_f = alpha_f * (1 / 0.7 - 1)

# Default-Priors for SIR-Model
DEFAULT_PRIORS = {
    "LAMBDA_0": {
        "type": "prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "params": {"mean": np.log(0.8), "sigma": 0.5},
        "name": r"$\lambda_0$",
        "description": "infectionrate before interventions"
    },
    "MU": {
        "type": "prior",
        "parameter_class": UnboundParameter,
        "distribution": np.random.lognormal,
        "params": {"mean": np.log(8), "sigma": 0.2},
        "name": r"$\mu$",
        "description": "recoveryrate"
    },
    "I0": {
        "type": "prior",
        "parameter_class": UnboundParameter,
        "distribution": np.random.gamma,
        # Expected value of gamma-distribution (shape * scale) = 60
        "params": {"shape": 2, "scale": 30},
        "name": r"$I_0$",
        "description": "initial number of infected"
    },
    "T1": {
        "type": "prior",
        "parameter_class": InterventionChangePointParameter,
        "distribution": np.random.normal,
        "params": {"loc": 8, "scale": 3},
        "name": r"$t_1$",
        "description": "first InterventionChangePoint"
    },
    "T2": {
        "type": "prior",
        "parameter_class": InterventionChangePointParameter,
        "distribution": np.random.normal,
        "params": {"loc": 15, "scale": 1},
        "name": r"$t_2$",
        "description": "second InterventionChangePoint"
    },
    "T3": {
        "type": "prior",
        "parameter_class": InterventionChangePointParameter,
        "distribution": np.random.normal,
        "params": {"loc": 22, "scale": 1},
        "name": r"$t_3$",
        "description": "third InterventionChangePoint"
    },
    "T4": {
        "type": "prior",
        "parameter_class": InterventionChangePointParameter,
        "distribution": np.random.normal,
        "params": {"loc": 66, "scale": 1},
        "name": r"$t_4$",
        "description": "fourth InterventionChangePoint"
    },
    "LAMBDA_1": {
        "type": "prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "params": {"mean": np.log(0.6), "sigma": 0.5},
        "name": r"$\lambda_1$",
        "description": "infectionrate after first intervention"
    },
    "LAMBDA_2": {
        "type": "prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "params": {"mean": np.log(0.3), "sigma": 0.5},
        "name": r"$\lambda_2$",
        "description": "infectionrate after second intervention"
    },
    "LAMBDA_3": {
        "type": "prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "params": {"mean": np.log(0.1), "sigma": 0.5},
        "name": r"$\lambda_3$",
        "description": "infectionrate after third intervention"
    },
    "LAMBDA_4": {
        "type": "prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "params": {"mean": np.log(0.1), "sigma": 0.5},
        "name": r"$\lambda_4$",
        "description": "infectionrate after fourth intervention"
    },
    "F_I": {
        "type": "prior",
        "parameter_class": WeeklyModulationAmplitudeParameter,
        "distribution": np.random.beta,
        # Erwartungswert der Beta-Verteilung alpha_f / (alpha_f + beta_f)
        "params": {"a": alpha_f, "b": beta_f},
        "name": r"$f_i$",
        "description": "modulation amplitude"
    },
    "PHI_I": {
        "type": "prior",
        "parameter_class": UnboundParameter,
        "distribution": stats.vonmises(kappa=0.01).rvs,
        "params": {},
        "name": r"$\phi_i$",
        "description": "phase shift"
    },
    "D_I": {
        "type": "prior",
        "parameter_class": DelayParameter,
        "distribution": np.random.lognormal,
        "params": {"mean": np.log(8), "sigma": 0.2},
        "name": r"$D_i$",
        "description": "reporting delay"
    },
    "PSI": {
        "type": "prior",
        "parameter_class": ScaleMultiplicativeReportingNoiseParameter,
        "distribution": np.random.gamma,
        "params": {"shape": 1, "scale": 5},
        "name": r"$\psi$",
        "description": "reporting noise"
    },
}


def create_prior_from_dd(key, additional_params={}):
    prior_entry = DEFAULT_PRIORS[key]
    return prior_entry["parameter_class"](distribution=partial(
        prior_entry["distribution"], **prior_entry["params"]), name=prior_entry["name"], **additional_params)


class SIRStrategy(ModelStrategy):
    @staticmethod
    def add_base(prior_array: list[UnboundParameter]) -> None:
        """
        Add base parameters to the prior array.
        """
        prior_array.append(create_prior_from_dd("LAMBDA_0"))
        prior_array.append(create_prior_from_dd("MU"))
        prior_array.append(create_prior_from_dd("I0"))

    @staticmethod
    def add_intervention(prior_array: list[UnboundParameter]) -> None:
        """
        Add intervention-related parameters to the prior array.
        """

        # could also use custom values from test_different_damping_priors
        prior_array.append(create_prior_from_dd("T1"))
        prior_array.append(create_prior_from_dd("T2"))
        prior_array.append(create_prior_from_dd("T3"))
        prior_array.append(create_prior_from_dd("T4"))
        prior_array.append(create_prior_from_dd("LAMBDA_1"))
        prior_array.append(create_prior_from_dd("LAMBDA_2"))
        prior_array.append(create_prior_from_dd("LAMBDA_3"))
        prior_array.append(create_prior_from_dd("LAMBDA_4"))

    @staticmethod
    def add_observation(prior_array: list[UnboundParameter], sim_lag: int) -> None:
        """
        Add observation-related parameters to the prior array.
        """
        prior_array.append(create_prior_from_dd("F_I"))
        prior_array.append(create_prior_from_dd("PHI_I"))
        prior_array.append(create_prior_from_dd(
            "D_I", additional_params={"sim_lag": sim_lag}))
        prior_array.append(create_prior_from_dd("PSI"))


def simulator_SIR(params: list[float], N: int, T: int, intervention_model: bool, observation_model: bool, param_names: list[str], fixed_params: dict[str, float] = dict(), sim_lag: int = 15) -> npt.NDArray[np.float64]:
    """Performs a forward simulation from the stationary SIR model given a random draw from the prior."""

    # Convert params array to a dictionary using param_names
    dynamic_params = {name: params[i] for i, name in enumerate(param_names)}

    # Merge with fixed parameters, overriding dynamic params where fixed
    combined_params = {**dynamic_params, **fixed_params}

    # Extract parameters using the Enum
    lambd0 = combined_params[DEFAULT_PRIORS["LAMBDA_0"]["name"]]
    mu = combined_params[DEFAULT_PRIORS["MU"]["name"]]
    I0 = combined_params[DEFAULT_PRIORS["I0"]["name"]]
    # Should IO also be a constraint instead of set to bouandry?
    I0 = max(1, np.round(I0))

    if intervention_model:
        t1 = combined_params[DEFAULT_PRIORS["T1"]["name"]]
        t2 = combined_params[DEFAULT_PRIORS["T2"]["name"]]
        t3 = combined_params[DEFAULT_PRIORS["T3"]["name"]]
        t4 = combined_params[DEFAULT_PRIORS["T4"]["name"]]
        # Round integer parameters
        # t1 = int(round(t1))
        t1, t2, t3, t4 = int(round(t1)), int(
            round(t2)), int(round(t3)), int(round(t4))

        lambd1 = combined_params[DEFAULT_PRIORS["LAMBDA_1"]["name"]]
        lambd2 = combined_params[DEFAULT_PRIORS["LAMBDA_2"]["name"]]
        lambd3 = combined_params[DEFAULT_PRIORS["LAMBDA_3"]["name"]]
        lambd4 = combined_params[DEFAULT_PRIORS["LAMBDA_4"]["name"]]

    if observation_model:
        f_i = combined_params[DEFAULT_PRIORS["F_I"]["name"]]
        phi_i = combined_params[DEFAULT_PRIORS["PHI_I"]["name"]]
        D_i = combined_params[DEFAULT_PRIORS["D_I"]["name"]]
        scale_I = combined_params[DEFAULT_PRIORS["PSI"]["name"]]
        D_i = int(round(D_i))

    # Maybe check for missing parameter

    # Configure model
    num_groups = 1
    sir_model = Model(num_groups)
    A0 = AgeGroup(0)

    # sim_lag is the maximum number of days for the delay
    # we simulat sim_lag days before t_0 to have values for t_0 - D
    t_max = T + sim_lag

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

    # Calculate damping value
    # Lambda0 is the initial contact rate which will be consecutively
    # reduced via the government measures
    def calc_damping_value_from_lambda(_lambdt: float, _lambd0: float = 1, _baseline: float = 1, _minimum: float = 0) -> float:
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

    if intervention_model:
        # delta_t1, delta_t2, delta_t3, delta_t4 = int(round(delta_t1)), int(round(delta_t2)), int(round(delta_t3)), int(round(delta_t4))
        sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
            coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd1), t=t1+sim_lag, level=0, type=0))
        sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
            coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd2), t=t2+sim_lag, level=0, type=0))
        sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
            coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd3), t=t3+sim_lag, level=0, type=0))
        sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(Damping(
            coeffs=np.ones((num_groups, num_groups)) * calc_damping_value_from_lambda(lambd4), t=t4+sim_lag, level=0, type=0))

    # Check logical constraints to parameters, should we really use apply_constraints?!
    sir_model.apply_constraints()

    # Reported new cases
    I_data = np.zeros(T)
    fs_i = np.zeros(T)

    # Run Simulation
    integrator = RKIntegratorCore(dt_max=1)
    (result, flows) = simulate_flows(
        0, t_max, 1, sir_model, integrator)

    # interpolate results
    flows = interpolate_simulation_result(flows)

    # consistency check, None values get deleted in configure input
    try:
        assert flows.get_last_time() == t_max
    except AssertionError as e:
        print('Invalid value simulated...return nan')
        return np.stack(([np.nan] * T, )).T

    if observation_model:

        # Adding new cases with delay D
        # Note, we assume the same delay
        shifted_t0 = sim_lag-D_i+1
        I_data = np.diff(flows.as_ndarray()[
                         1, shifted_t0:shifted_t0+T+1])
        I_data = np.clip(I_data, 10 ** -14, N)

        # Compute lags
        fs_i = (1-f_i)*(1 -
                        np.abs(np.sin((np.pi/7) * np.arange(0, T, 1) - 0.5*phi_i)))

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
        DEFAULT_PRIORS["LAMBDA_0"]["name"], DEFAULT_PRIORS["MU"]["name"], DEFAULT_PRIORS["I0"]["name"]])

    plt.plot(results)
    plt.savefig("simulator_SIR.png")
