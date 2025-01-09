import numpy as np
import numpy.typing as npt
from scipy import stats
from functools import partial
import pandas as pd
import os
from copy import deepcopy

import memilio.simulation as mio
import memilio.simulation.osir as osir
from memilio.inference.prior import (UnboundParameter, LambdaParameter, DelayParameter,
                                     InterventionChangePointParameter, WeeklyModulationAmplitudeParameter,
                                     ScaleMultiplicativeReportingNoiseParameter, ModelStrategy, BoundParameter)
from memilio.inference.networks import NETWORKS_DICT

from bayesflow.simulation import TwoLevelPrior, ContextGenerator
import matplotlib.pyplot as plt

alpha_f = (0.7**2)*((1-0.7)/(0.17**2) - (1-0.7))
beta_f = alpha_f*(1/0.7 - 1)
mio.set_log_level(mio.LogLevel.Warning)

# Default-Priors for SIR-MLM-Model
DEFAULT_PRIORS = {
    "LAMBDA_0_MEAN": {
        "type": "hyper_prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.normal,
        "dist_params": {"loc": 0.8, "scale": 0.5},
        "name": r"$\lambda_0_mean$",
        "description": "mean of lambda_0 prior"
    },
    "LAMBDA_0_SIGMA": {
        "type": "hyper_prior",
        "parameter_class": BoundParameter,
        "parameter_class_params": {"lower_bound": 0},
        "distribution": np.random.normal,
        "dist_params": {"loc": 0.5, "scale": 0.2},
        "name": r"$\lambda_0_sigma$",
        "description": "sigma of lambda_0 prior"
    },
    "LAMBDA_0": {
        "type": "local_prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "dist_params": {},
        "name": r"$\lambda_0$",
        "description": "infectionrate before interventions"
    },
    "I0": {
        "type": "local_prior",
        "parameter_class": UnboundParameter,
        "distribution": np.random.gamma,
        # Expected value of gamma-distribution (shape * scale) = 60
        "dist_params": {"shape": 2, "scale": 30},
        "name": r"$I_0$",
        "description": "initial number of infected"
    },
    "MU": {
        "type": "shared_prior",
        "parameter_class": UnboundParameter,
        "distribution": np.random.lognormal,
        "dist_params": {"mean": np.log(8), "sigma": 0.2},
        "name": r"$\mu$",
        "description": "recoveryrate"
    },
    "T1": {
        "type": "prior",
        "parameter_class": InterventionChangePointParameter,
        "distribution": np.random.normal,
        "dist_params": {"loc": 8, "scale": 3},
        "name": r"$t_1$",
        "description": "first InterventionChangePoint"
    },
    "T2": {
        "type": "prior",
        "parameter_class": InterventionChangePointParameter,
        "distribution": np.random.normal,
        "dist_params": {"loc": 15, "scale": 1},
        "name": r"$t_2$",
        "description": "second InterventionChangePoint"
    },
    "T3": {
        "type": "prior",
        "parameter_class": InterventionChangePointParameter,
        "distribution": np.random.normal,
        "dist_params": {"loc": 22, "scale": 1},
        "name": r"$t_3$",
        "description": "third InterventionChangePoint"
    },
    "T4": {
        "type": "prior",
        "parameter_class": InterventionChangePointParameter,
        "distribution": np.random.normal,
        "dist_params": {"loc": 66, "scale": 1},
        "name": r"$t_4$",
        "description": "fourth InterventionChangePoint"
    },
    "LAMBDA_1": {
        "type": "prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "dist_params": {"mean": np.log(0.6), "sigma": 0.5},
        "name": r"$\lambda_1$",
        "description": "infectionrate after first intervention"
    },
    "LAMBDA_2": {
        "type": "prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "dist_params": {"mean": np.log(0.3), "sigma": 0.5},
        "name": r"$\lambda_2$",
        "description": "infectionrate after second intervention"
    },
    "LAMBDA_3": {
        "type": "prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "dist_params": {"mean": np.log(0.1), "sigma": 0.5},
        "name": r"$\lambda_3$",
        "description": "infectionrate after third intervention"
    },
    "LAMBDA_4": {
        "type": "prior",
        "parameter_class": LambdaParameter,
        "distribution": np.random.lognormal,
        "dist_params": {"mean": np.log(0.1), "sigma": 0.5},
        "name": r"$\lambda_4$",
        "description": "infectionrate after fourth intervention"
    },
    "F_I": {
        "type": "shared_prior",
        "parameter_class": WeeklyModulationAmplitudeParameter,
        "distribution": np.random.beta,
        # Erwartungswert der Beta-Verteilung alpha_f / (alpha_f + beta_f)
        "dist_params": {"a": alpha_f, "b": beta_f},
        "name": r"$f_i$",
        "description": "modulation amplitude"
    },
    "PHI_I": {
        "type": "shared_prior",
        "parameter_class": UnboundParameter,
        "distribution": stats.vonmises(kappa=0.01).rvs,
        "dist_params": {},
        "name": r"$\phi_i$",
        "description": "phase shift"
    },
    "D_I": {
        "type": "shared_prior",
        "parameter_class": DelayParameter,
        "distribution": np.random.lognormal,
        "dist_params": {"mean": np.log(8), "sigma": 0.2},
        "name": r"$D_i$",
        "description": "reporting delay"
    },
    "PSI": {
        "type": "shared_prior",
        "parameter_class": ScaleMultiplicativeReportingNoiseParameter,
        "distribution": np.random.gamma,
        "dist_params": {"shape": 1, "scale": 5},
        "name": r"$\psi$",
        "description": "reporting noise"
    },
}


def get_param_names_by_type(default_priors, param_type):
    """! Retrieves all parameter names from the DEFAULT_PRIORS dictionary based on their type.

    @param param_type The type of parameters to filter as str (e.g., 'local_prior', 'shared_prior').
    """
    return [key for key, value in default_priors.items() if value.get("type") == param_type]


def create_prior_from_dd(key, additional_params={}, hyper_params={}):
    prior_entry = DEFAULT_PRIORS[key]
    return prior_entry["parameter_class"](distribution=partial(
        prior_entry["distribution"], **prior_entry.get("dist_params", {}), **hyper_params), name=prior_entry["name"], **prior_entry.get("parameter_class_params", {}), **additional_params)


def create_prior():

    def draw_hyper_prior():
        # Draw location for 2D conditional prior
        lambd0_mean = create_prior_from_dd("LAMBDA_0_MEAN").get_draw()()
        lambd0_mean = np.log(lambd0_mean)  # np.log(0.8)
        lambd0_sigma = create_prior_from_dd(
            "LAMBDA_0_SIGMA").get_draw()()  # 0.5
        return [lambd0_mean, lambd0_sigma]

    def draw_local_prior(hyper_params, n_regions):
        lambd0_mean, lambd0_sigma = hyper_params

        lambd0 = np.zeros(n_regions)
        I0 = np.zeros(n_regions)
        for i in range(n_regions):
            # Draw parameter given location from hyperprior
            lambd0[i] = create_prior_from_dd("LAMBDA_0", hyper_params={
                                             "mean": lambd0_mean, "sigma": lambd0_sigma}).get_draw()()
            # every state has its own I0
            I0[i] = create_prior_from_dd("I0").get_draw()()

        return [lambd0, I0]

    def draw_shared_prior(sim_lag=15):

        mu = create_prior_from_dd("MU").get_draw()()

        # observation model
        f_i = create_prior_from_dd("F_I").get_draw()()
        phi_i = create_prior_from_dd("PHI_I").get_draw()()
        d_i = create_prior_from_dd(
            "D_I", additional_params={"sim_lag": sim_lag}).get_draw()()
        scale_I = create_prior_from_dd("PSI").get_draw()()
        return [mu, f_i, phi_i, d_i, scale_I]

    # context = ContextGenerator(non_batchable_context_fun=lambda : np.random.randint(1, 101))
    # prior = TwoLevelPrior(draw_hyper_prior, draw_local_prior, draw_shared_prior, local_context_generator = context)

    # prior(batch_size = 320, local_args = {n_regions: n_groups})

    return TwoLevelPrior(
        draw_hyper_prior, draw_local_prior, draw_shared_prior)


def build_mlm_amortizer(_settings: dict, local_params, global_params):
    # # Should not be lower than number of parameters
    # summary_net = SequenceNetwork(summary_dim=trainer_parameters.summary_dim)

    model_settings = deepcopy(_settings["model"])

    # could try setting coupling_design=spline for more superior performance on lower-dimensional problems.
    # infere local parameters on local data
    local_inference_net = NETWORKS_DICT[model_settings["local_amortizer"]["inference_net"].pop(
        "type")](num_params=local_params, **model_settings["local_amortizer"]["inference_net"])
    local_amortizer = NETWORKS_DICT[model_settings["local_amortizer"]["type"]](
        local_inference_net, name=model_settings["local_amortizer"]["name"])

    # could try setting coupling_design=spline for more superior performance on lower-dimensional problems.
    # infere hyper and shared parameters on global data
    global_inference_net = NETWORKS_DICT[model_settings["global_amortizer"]["inference_net"].pop(
        "type")](num_params=global_params, **model_settings["global_amortizer"]["inference_net"])
    global_amortizer = NETWORKS_DICT[model_settings["global_amortizer"]["type"]](
        global_inference_net, name=model_settings["global_amortizer"]["name"])

    local_summary_net = NETWORKS_DICT[model_settings["summary_net_kwargs"]["local_summary_net"].pop(
        "type")](**model_settings["summary_net_kwargs"]["local_summary_net"])
    global_summary_net = NETWORKS_DICT[model_settings["summary_net_kwargs"]["global_summary_net"].pop(
        "type")](**model_settings["summary_net_kwargs"]["global_summary_net"])
    summary_net = NETWORKS_DICT[model_settings["summary_net_kwargs"]["type"]](
        [local_summary_net, global_summary_net])

    amortizer = NETWORKS_DICT[model_settings["type"]](
        local_amortizer, global_amortizer, summary_net, name=model_settings["name"])

    return amortizer


State = {
    1: 'Schleswig-Holstein',
    2: 'Hamburg',
    3: 'Niedersachsen',
    4: 'Bremen',
    5: 'Nordrhein-Westfalen',
    6: 'Hessen',
    7: 'Rheinland-Pfalz',
    8: 'Baden-Württemberg',
    9: 'Bayern',
    10: 'Saarland',
    11: 'Berlin',
    12: 'Brandenburg',
    13: 'Mecklenburg-Vorpommern',
    14: 'Sachsen',
    15: 'Sachsen-Anhalt',
    16: 'Thüringen',


}


def get_population(path):
    """! read population data in list from dataset
    @param path Path to the dataset containing the population data
    """

    data = pd.read_json(path)
    data["ID_State"] = data["ID_County"] // 1000
    population = data.groupby("ID_State")["Population"].sum()

    return np.array(population)


def simulator_SIR(combined_params: tuple[list[list[float]], list[float]], n_regions: int, population: list[int], T: int, mobility_params: list[list[int]], intervention_model: bool, observation_model: bool, param_names: list[str], sim_lag: int = 15) -> npt.NDArray[np.float64]:
    """Performs a forward simulation from the stationary SIR model given a random draw from the prior."""

    (local_params, shared_params) = combined_params
    lambd0, I0 = local_params
    mu, f_i, phi_i, D_i, scale_I = shared_params
    I0 = np.fmax(20, np.round(I0))
    D_i = int(round(D_i))

    # if intervention_model:
    #     t1 = params[DEFAULT_PRIORS["T1"]["name"]]
    #     t2 = params[DEFAULT_PRIORS["T2"]["name"]]
    #     t3 = params[DEFAULT_PRIORS["T3"]["name"]]
    #     t4 = params[DEFAULT_PRIORS["T4"]["name"]]
    #     # Round integer parameters
    #     t1, t2, t3, t4 = int(round(t1)), int(
    #         round(t2)), int(round(t3)), int(round(t4))

    #     lambd1 = params[DEFAULT_PRIORS["LAMBDA_1"]["name"]]
    #     lambd2 = params[DEFAULT_PRIORS["LAMBDA_2"]["name"]]
    #     lambd3 = params[DEFAULT_PRIORS["LAMBDA_3"]["name"]]
    #     lambd4 = params[DEFAULT_PRIORS["LAMBDA_4"]["name"]]

    # if observation_model:
    #     f_i = params[DEFAULT_PRIORS["F_I"]["name"]]
    #     phi_i = params[DEFAULT_PRIORS["PHI_I"]["name"]]
    #     D_i = params[DEFAULT_PRIORS["D_I"]["name"]]
    #     scale_I = params[DEFAULT_PRIORS["PSI"]["name"]]
    #     D_i = int(round(D_i))

    # Maybe check for missing parameter

    # check mobility_params dim and 0 values at diagonal

    # sim_lag is the maximum number of days for the delay
    # we simulat sim_lag days before t_0 to have values for t_0 - D
    t_max = T + sim_lag

    # Configure model
    num_agegroups = 1
    graph = osir.MobilityGraph()

    for i in range(n_regions):
        sir_model = osir.Model(num_agegroups=num_agegroups)
        A0 = mio.AgeGroup(0)

        # Initial conditions
        sir_model.populations[A0, osir.InfectionState.Infected] = I0[i]
        sir_model.populations[A0, osir.InfectionState.Recovered] = 0
        sir_model.populations.set_difference_from_total(
            (A0, osir.InfectionState.Susceptible), population[i])

        # Initialize Parameters
        sir_model.parameters.TimeInfected[A0] = mu
        sir_model.parameters.TransmissionProbabilityOnContact[A0] = 1
        # Very weird way of setting dampings and differnt from paper, because of no time till NPIs fully take effect
        # Also damping could increase the beta value, should that be possible?

        def calc_damping_value_from_lambda(_lambdt: float, _lambd0: float = 1, _baseline: float = 1, _minimum: float = 0) -> float:
            # cf = bl - (dampingvalue * (bl - min))
            # -> bl - cf = (dampingvalue * (bl - min)
            # -> (bl - cf) / (bl - min) = dampingvalue
            # lambdt = cf * lambd0 -> cf = lambdt / lambd0
            # => dampingvalue = (bl - (lambdt / lambd0)) / (bl - min)
            return (_baseline - (_lambdt / _lambd0)) / (_baseline - _minimum)

        sir_model.parameters.ContactPatterns.cont_freq_mat[0].baseline = np.ones(
            (num_agegroups, num_agegroups))
        sir_model.parameters.ContactPatterns.cont_freq_mat[0].minimum = np.zeros(
            (num_agegroups, num_agegroups))
        sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(mio.Damping(
            coeffs=np.ones((num_agegroups, num_agegroups)) * calc_damping_value_from_lambda(lambd0[i]), t=0, level=0, type=0))

        if intervention_model:
            # delta_t1, delta_t2, delta_t3, delta_t4 = int(round(delta_t1)), int(round(delta_t2)), int(round(delta_t3)), int(round(delta_t4))
            sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(mio.Damping(
                coeffs=np.ones((num_agegroups, num_agegroups)) * calc_damping_value_from_lambda(lambd1[i]), t=t1+sim_lag, level=0, type=0))
            sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(mio.Damping(
                coeffs=np.ones((num_agegroups, num_agegroups)) * calc_damping_value_from_lambda(lambd2[i]), t=t2+sim_lag, level=0, type=0))
            sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(mio.Damping(
                coeffs=np.ones((num_agegroups, num_agegroups)) * calc_damping_value_from_lambda(lambd3[i]), t=t3+sim_lag, level=0, type=0))
            sir_model.parameters.ContactPatterns.cont_freq_mat.add_damping(mio.Damping(
                coeffs=np.ones((num_agegroups, num_agegroups)) * calc_damping_value_from_lambda(lambd4[i]), t=t4+sim_lag, level=0, type=0))

        # Check logical constraints to parameters, should we really use apply_constraints?!
        sir_model.apply_constraints()

        graph.add_node(id=0, model=sir_model, t0=0)

        for end_node_idx in range(n_regions):
            mobility_coefficients = (mobility_params[i][end_node_idx] / population[i]) * \
                np.ones(sir_model.populations.numel())

            graph.add_edge(i, end_node_idx, mio.MobilityParameters(
                mobility_coefficients))

    # run simulation
    sim = osir.MobilitySimulation(graph, 0, dt=0.5)
    sim.advance(t_max)

    result = [osir.interpolate_simulation_result(sim.graph.get_node(
        i).property.result) for i in range(n_regions)]

    out_data = []

    if observation_model:
        fs_i = np.zeros(T)
        # Compute lags
        fs_i = (1-f_i)*(1 -
                        np.abs(np.sin((np.pi/7) * np.arange(0, T, 1) - 0.5*phi_i)))

        for i in range(n_regions):
            # Reported new cases
            I_data = np.zeros(T)

            # Adding new cases with delay D
            # Note, we assume the same delay
            I_data = np.diff(result[i].as_ndarray()[
                             1, (sim_lag-D_i):(sim_lag-D_i)+T+1]) * (-1)
            I_data = np.clip(I_data, 10 ** -14, population[i])

            # Compute weekly modulation
            I_data = (1-fs_i) * I_data

            # check for negative values
            try:
                scale = np.sqrt(I_data)*scale_I
                assert np.all(scale >= 0)
            except AssertionError as e:
                print('Invalid value simulated...return nan')
                return np.stack(([[np.nan] * T] * n_regions, )).T

            # Add noise
            I_data = stats.t(df=4, loc=I_data,
                             scale=np.sqrt(I_data)*scale_I).rvs()

            # bound all negative values to 0
            I_data = np.clip(I_data, 10 ** -14, population[i])

            out_data.append(I_data)
    else:
        for i in range(n_regions):
            # Reported new cases
            I_data = np.zeros(T)
            I_data = np.diff(result[i].as_ndarray()[1, :T+1]) * (-1)

            # bound all negative values to 0
            I_data = np.clip(I_data, 10 ** -14, population[i])

            out_data.append(I_data)

    return np.stack((out_data, )).transpose((1, 2, 0))

    # at this point code example uses cum sum instead of new
    # return np.stack((flows[:,0], )).T # np.stack((result[:,1])).T also wrong because its not cumsum


if __name__ == "__main__":

    path_population_data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                                        "..", "..", "data", "pydata", "Germany",
                                        "county_current_population.json")

    population = get_population(path_population_data)
    prior = create_prior()
    params = prior(1, local_args={"n_regions": 16})
    results = simulator_SIR([params["local_parameters"][0], params["shared_parameters"][0]], 16, population, 81, False, False, [
        DEFAULT_PRIORS["LAMBDA_0"]["name"], DEFAULT_PRIORS["MU"]["name"], DEFAULT_PRIORS["I0"]["name"]])

    plt.plot(results)
    plt.savefig("simulator_SIR.png")
