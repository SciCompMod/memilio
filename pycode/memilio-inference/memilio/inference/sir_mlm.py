import numpy as np
import numpy.typing as npt
from scipy import stats
from enum import Enum
from functools import partial
import json
import pandas as pd
import os

import memilio.simulation as mio
import memilio.simulation.osir as osir
from memilio.inference.prior import (UnboundParameter, LambdaParameter, DelayParameter,
                                     InterventionChangePointParameter, WeeklyModulationAmplitudeParameter,
                                     ScaleMultiplicativeReportingNoiseParameter, ModelStrategy, BoundParameter)

from bayesflow.simulation import TwoLevelPrior, ContextGenerator
import matplotlib.pyplot as plt

alpha_f = (0.7**2)*((1-0.7)/(0.17**2) - (1-0.7))
beta_f = alpha_f*(1/0.7 - 1)
mio.set_log_level(mio.LogLevel.Warning)


class HyperparameterNamesSir(Enum):
    MEAN = r'$mean$'
    SIGMA = r'$sigma$'


class ParameterNamesSir(Enum):
    LAMBDA_0 = r'$\lambda_0$'
    MU = r'$\mu$'
    I0 = r'$I_0$'
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
    D_I = r'$D_i$'
    PSI = r'$\psi$'


class SIRStrategy(ModelStrategy):
    @staticmethod
    # Possible option to draw without redraw
    def add_base(prior_array: list[UnboundParameter]) -> None:
        prior_array.append(LambdaParameter(distribution=partial(
            np.random.lognormal, mean=np.log(1.2), sigma=0.5), name=ParameterNamesSir.LAMBDA_0.value))
        prior_array.append(UnboundParameter(distribution=partial(
            np.random.lognormal, mean=np.log(8), sigma=0.2), name=ParameterNamesSir.MU.value))
        prior_array.append(UnboundParameter(distribution=partial(
            np.random.gamma, shape=2, scale=30), name=ParameterNamesSir.I0.value))

    @staticmethod
    # Possible option to draw without redraw
    def add_intervention(prior_array: list[UnboundParameter]) -> None:
        prior_array.append(InterventionChangePointParameter(distribution=partial(
            np.random.normal, loc=8, scale=3), name=ParameterNamesSir.T1.value))
        prior_array.append(InterventionChangePointParameter(distribution=partial(
            np.random.normal, loc=15, scale=1), name=ParameterNamesSir.T2.value))
        prior_array.append(InterventionChangePointParameter(distribution=partial(
            np.random.normal, loc=22, scale=1), name=ParameterNamesSir.T3.value))
        prior_array.append(InterventionChangePointParameter(distribution=partial(
            np.random.normal, loc=66, scale=1), name=ParameterNamesSir.T4.value))
        prior_array.append(LambdaParameter(distribution=partial(
            np.random.lognormal, mean=np.log(0.6), sigma=0.5), name=ParameterNamesSir.LAMBDA_1.value))
        prior_array.append(LambdaParameter(distribution=partial(
            np.random.lognormal, mean=np.log(0.3), sigma=0.5), name=ParameterNamesSir.LAMBDA_2.value))
        prior_array.append(LambdaParameter(distribution=partial(
            np.random.lognormal, mean=np.log(0.1), sigma=0.5), name=ParameterNamesSir.LAMBDA_3.value))
        prior_array.append(LambdaParameter(distribution=partial(
            np.random.lognormal, mean=np.log(0.1), sigma=0.5), name=ParameterNamesSir.LAMBDA_4.value))

    @staticmethod
    # Possible option to draw without redraw
    def add_observation(prior_array: list[UnboundParameter], sim_diff: int) -> None:
        prior_array.append(WeeklyModulationAmplitudeParameter(distribution=partial(
            np.random.beta, a=alpha_f, b=beta_f), name=ParameterNamesSir.F_I.value))
        prior_array.append(UnboundParameter(
            distribution=stats.vonmises(kappa=0.01).rvs, name=ParameterNamesSir.PHI_I.value))
        prior_array.append(DelayParameter(distribution=partial(
            np.random.lognormal, mean=np.log(8), sigma=0.2), name=ParameterNamesSir.D_I.value, sim_diff=sim_diff))
        prior_array.append(ScaleMultiplicativeReportingNoiseParameter(distribution=partial(
            np.random.gamma, shape=1, scale=5), name=ParameterNamesSir.PSI.value))


def create_prior():

    def draw_hyper_prior():
        # Draw location for 2D conditional prior
        mean = LambdaParameter(distribution=partial(
            np.random.normal, loc=0.8, scale=0.2), name=HyperparameterNamesSir.MEAN.value).get_draw()()  # np.log(0.8)
        sigma = BoundParameter(distribution=partial(
            np.random.normal, loc=0.5, scale=0.2), name=HyperparameterNamesSir.SIGMA.value, lower_bound=0).get_draw()()  # 0.5
        return [mean, sigma]

    def draw_local_prior(hyper_params, n_regions):
        mean, sigma = hyper_params

        lambd0 = np.zeros(n_regions)
        I0 = np.zeros(n_regions)
        for i in range(n_regions):
            # Draw parameter given location from hyperprior
            lambd0[i] = LambdaParameter(distribution=partial(
                np.random.lognormal, mean=np.log(mean), sigma=sigma), name=ParameterNamesSir.LAMBDA_0.value).get_draw()()

            # every state has its own I0
            I0[i] = UnboundParameter(distribution=partial(
                np.random.gamma, shape=2, scale=30), name=ParameterNamesSir.I0.value).get_draw()()

        return [lambd0, I0]

    def draw_shared_prior(sim_lag=15):

        mu = UnboundParameter(distribution=partial(
            np.random.lognormal, mean=np.log(8), sigma=0.2), name=ParameterNamesSir.MU.value).get_draw()()

        # observation model
        f_i = WeeklyModulationAmplitudeParameter(distribution=partial(
            np.random.beta, a=alpha_f, b=beta_f), name=ParameterNamesSir.F_I.value).get_draw()()
        phi_i = UnboundParameter(
            distribution=stats.vonmises(kappa=0.01).rvs, name=ParameterNamesSir.PHI_I.value).get_draw()()
        d_i = DelayParameter(distribution=partial(np.random.lognormal, mean=np.log(
            8), sigma=0.2), name=ParameterNamesSir.D_I.value, sim_lag=sim_lag).get_draw()()
        scale_I = ScaleMultiplicativeReportingNoiseParameter(distribution=partial(
            np.random.gamma, shape=1, scale=5), name=ParameterNamesSir.PSI.value).get_draw()()
        return [mu, f_i, phi_i, d_i, scale_I]

    # context = ContextGenerator(non_batchable_context_fun=lambda : np.random.randint(1, 101))
    # prior = TwoLevelPrior(draw_hyper_prior, draw_local_prior, draw_shared_prior, local_context_generator = context)

    # prior(batch_size = 320, local_args = {n_regions: n_groups})

    return TwoLevelPrior(
        draw_hyper_prior, draw_local_prior, draw_shared_prior)


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


def simulator_SIR(combined_params: tuple[list[list[float]], list[float]], n_regions: int, N: list[int], T: int, mobility_params: list[list[int]], intervention_model: bool, observation_model: bool, param_names: list[str], sim_lag: int = 15) -> npt.NDArray[np.float64]:
    """Performs a forward simulation from the stationary SIR model given a random draw from the prior."""

    (local_params, shared_params) = combined_params
    lambd0, I0 = local_params
    mu, f_i, phi_i, D_i, scale_I = shared_params
    I0 = np.fmax(20, np.round(I0))
    D_i = int(round(D_i))

    # if intervention_model:
    #     t1 = params[ParameterNamesSir.T1.value]
    #     t2 = params[ParameterNamesSir.T2.value]
    #     t3 = params[ParameterNamesSir.T3.value]
    #     t4 = params[ParameterNamesSir.T4.value]
    #     # Round integer parameters
    #     t1, t2, t3, t4 = int(round(t1)), int(
    #         round(t2)), int(round(t3)), int(round(t4))

    #     lambd1 = params[ParameterNamesSir.LAMBDA_1.value]
    #     lambd2 = params[ParameterNamesSir.LAMBDA_2.value]
    #     lambd3 = params[ParameterNamesSir.LAMBDA_3.value]
    #     lambd4 = params[ParameterNamesSir.LAMBDA_4.value]

    # if observation_model:
    #     f_i = params[ParameterNamesSir.F_I.value]
    #     phi_i = params[ParameterNamesSir.PHI_I.value]
    #     D_i = params[ParameterNamesSir.D_I.value]
    #     scale_I = params[ParameterNamesSir.PSI.value]
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
            (A0, osir.InfectionState.Susceptible), N[i])

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
            mobility_coefficients = (mobility_params[i][end_node_idx] / N[i]) * \
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
            I_data = np.clip(I_data, 10 ** -14, N[i])

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
            I_data = np.clip(I_data, 10 ** -14, N[i])

            out_data.append(I_data)
    else:
        for i in range(n_regions):
            # Reported new cases
            I_data = np.zeros(T)
            I_data = np.diff(result[i].as_ndarray()[1, :T+1]) * (-1)

            # bound all negative values to 0
            I_data = np.clip(I_data, 10 ** -14, N[i])

            out_data.append(I_data)

    return np.stack((out_data, )).transpose((1, 2, 0))

    # at this point code example uses cum sum instead of new
    # return np.stack((flows[:,0], )).T # np.stack((result[:,1])).T also wrong because its not cumsum


if __name__ == "__main__":

    path_population_data = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..",
                                        "..", "..", "data", "pydata", "Germany",
                                        "county_current_population.json")

    N = get_population(path_population_data)
    prior = create_prior()
    params = prior(1, local_args={"n_regions": 16})
    results = simulator_SIR([params["local_parameters"][0], params["shared_parameters"][0]], 16, N, 81, False, False, [
        ParameterNamesSir.LAMBDA_0.value, ParameterNamesSir.MU.value, ParameterNamesSir.I0.value])

    plt.plot(results)
    plt.savefig("simulator_SIR.png")
